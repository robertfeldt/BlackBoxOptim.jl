log_k_temp(k, T0) = T0 / log(k)

gaussian_sample_neighbor(s, T, D, ranges) = s + T .* randn(D)

function asa_sample_neighbor(s, T, D, ranges)
  u = rand(D)
  y = sign(u - 0.5) .* T .* ((1 + 1 ./ T).^abs(2*u-1) - 1)
  s .+ (y .* ranges)
end

function simulated_annealing(func,
                             s0,
                             numiterations = 10000,
                             sample_neighbor = asa_sample_neighbor,
                             temperature = log_k_temp,
                             ranges = ones(length(s0)))
 
  sbest = s = s0
  es = func(s)
  ebest = func(sbest)
  T0 = 100.0
  D = length(s)

  for k = 1:numiterations

    t = log_k_temp(k, T0)

    snew = sample_neighbor(s, t, D, ranges)
    enew = func(snew)

    delta_e = enew - es

    if delta_e < 0

      s = snew
      es = enew
      if es < ebest
        ebest = es
        sbest = s
        println("$(k), New best: $(ebest)")
      end

    else

      if rand() <= exp( -delta_e / t )
        s = snew
        es = enew
      end

    end

  end

  return sbest, ebest

end

function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

simulated_annealing(rosenbrock2d, randn(2), 100000)