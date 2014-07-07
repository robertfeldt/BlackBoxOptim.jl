@everywhere function fake_slow_condition(x, delay = 0.10, expected = 30)
  sleep(delay)
  sum(x) == expected ? x : false # Return x iff it fulfills the poperty
end

# Similar to Base.findfirst(predicate, A) but process the array "A" columnwise in parallel
# and opt out of further processing if found. Returns both the index at which found
# as well as the result returned.
@everywhere function pfind_colwise(predicate, A)
  np = nprocs()  # determine the number of processes available
  n = size(A, 2)
  i = 1
  results = Any[false for i in 1:n]
  found = false
  wasfound(idx) = found = true
  isfound() = found
  # function to produce the next work item from the queue.
  # in this case it's just an index to a column.
  nextidx() = (idx=i; i+=1; idx)
  @sync begin
    for p=1:np
      if p != myid() || np == 1
        @async begin
          while true
            idx = nextidx()
            if idx > n
              break
            end
            results[idx] = remotecall_fetch(p, predicate, A[:,idx])
            if isfound()
              break # Someone else already found it
            end
            if results[idx] != false
              wasfound(idx)
              break # We found it
            end
          end
        end
      end
    end
  end
  fidx = findfirst((r)->(r != false), results)
  (fidx, (fidx != 0 ? results[fidx] : nothing))
end
