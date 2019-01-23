# This is a general type for Adaptive selection based on Area-under-curve
# multi-armed bandit as described in the paper:
#   Francois-Michel De Rainville et al, "Sustainable Cooperative Coevolution 
#   with a Multi-Armed Bandit", 2013
#

# An AdaptiveSelector selects between a set of integers that represent different
# operators, values etc. It does not have any info about what the things it selects
# between are, below we call them items. The set of items the selector selects
# from is typically static but might not be.
abstract type AdaptiveSelector end

# The default AdaptiveSelector is a random selector.
mutable struct RandomSelector <: AdaptiveSelector
  item_map::Dict{Int, Any}
  next_index::Int

  RandomSelector(items = []) = begin
    rs = new(Dict{Int, Any}(), 0)
    add_items(rs, items)
    rs
  end
end

# Returns an item map that maps from the indices of the current item set to
# each item.
function item_map(a::AdaptiveSelector)
  a.item_map # Override this if you do not save the item map in an instance var
end

# Get next free index. Returns unique index number each time it is called.
function get_next_free_index(a::AdaptiveSelector)
  a.next_index += 1
end

# Add an item. Returns its index.
function add_item(a::AdaptiveSelector, item)
  index = get_next_free_index(a)
  item_map(a)[index] = item
  index
end

function add_items(a::AdaptiveSelector, items)
  [add_item(a, item) for item in items]
end

# Remove an item given its index.
function remove_item_with_index(a::AdaptiveSelector, item_index::Int)
  delete!(a.item_map, item_index)
end

# Remove an item.
function remove_item(a::AdaptiveSelector, item)
  not_equals_item = (k,v) -> v != item
  filter!( not_equals_item, a.item_map )
end


# Select one item in the set of things that we select from. Returns both its
# index and itself / its value.
function select(a::AdaptiveSelector)
  ks = keys(a.item_map)
  randkey = ks[rand(1:length(ks))]
  a.item_map[randkey]
end

# Add reward for one item given its index. Default is to not do anything since
# the default selector is a random selector that need not care about rewards.
# A reward should be a positive value with larger values indicating more 
# benefit when using the item.
function reward(a::AdaptiveSelector, index, reward)
  # Override for selectors that use rewards.
end

# We need a queue of fixed size to keep info about recent rewards.
mutable struct FixedSizeQueue
  queue::Array{Any,1}
  size::Int
  maxsize::Int
  index::Int

  FixedSizeQueue(maxsize) = begin
    q = Any[]
    [Base.push!(q, 0) for i in 1:maxsize]
    new(q, 0, maxsize, 0)
  end
end

size(fsq::FixedSizeQueue) = fsq.size
maxsize(fsq::FixedSizeQueue) = fsq.maxsize
is_full(fsq::FixedSizeQueue) = fsq.size == fsq.maxsize

function push!(fsq::FixedSizeQueue, item)
  fsq.queue[fsq.index + 1] = item
  fsq.index = mod(fsq.index + 1, fsq.maxsize)
  fsq.size = min(fsq.maxsize, fsq.size + 1)
end

function as_vector(fsq::FixedSizeQueue)
  v = Any[]
  for elem = fsq
    push!(v, elem)
  end
  v
end

##############
# Define FixedSizeQueue to be an iterator by defining start, done and next:
##############

Base.start(fsq::FixedSizeQueue) = (fsq.size, mod(fsq.index - fsq.size, fsq.maxsize))

Base.done(fsq::FixedSizeQueue, state) = state[1] == 0

function Base.next(fsq::FixedSizeQueue, state)
  elem = fsq.queue[state[2] + 1]
  return elem, (state[1]-1, mod(state[2]+1, fsq.maxsize))
end


mutable struct AucBanditSelector <: AdaptiveSelector
  window::FixedSizeQueue  # Window of latest (itemindex, reward) pairs
  ns::Int[]             # Number of times each item has been used within the current window.
  items::Any[]
  next_index::Int
  decay_factor::Float64
  exploration_factor::Float64

  AucBanditSelector(items = []; window_size = 50, decay_factor = 1.0,
    exploration_factor = 1.0) = begin

    aucb = new(FixedSizeQueue(window_size), Int[], Any[], 0, 
      decay_factor, exploration_factor)

    add_items(aucb, items)

    aucb

  end
end

function add_item(a::AucBanditSelector, item)
  push!(a.ns, 0) # Add a zero to indicate how often it has been used => none!
  push!(a.items, item)
  length(a.items)
end

# Remove an item given its index.
function remove_item_with_index(a::AucBanditSelector, item_index::Int)
  splice!(a.items, item_index)
  splice!(a.ns, item_index)
end

# Remove an item.
function remove_item(a::AucBanditSelector, item)
  for i in 1:length(a.items)
    if a.items[i] == item
      remove_item_with_index(a, i)
      return i
    end
  end
end

function reward(a::AucBanditSelector, index, reward)
  # Update the ns counts by subtracting one from the one that is about the exit 
  # the window.
  if is_full(a.window)
    (index, reward) = first(a.window)
    a.ns[index] -= 1
  end
  push!(a.window, (index, reward))
  a.ns[index] += 1
end

function select(a::AucBanditSelector)
  unused = a.ns == 0
  if any(unusued)
    # Randomly select one of the unused
    indices_to_unused = collect(1:length(a.ns))[unused]
    indices_to_unused[rand(1:length(indices_to_unused))]
  else 
    # Select via the multi-armed bandit formula
    qs = calc_quality_values(a)
    vs = qs + a.exploration_factor * sqrt( 2 * log(sum(a.ns)) * (1 / a.ns) )
    items[argmin(vs)]
  end
end

function calc_quality_values(a:AucBanditSelector)
  qs = zeros(length(a.ns))
  for i in 1:length(a.ns)
    qs[i] = area_under_curve_credit_assignment(a, i)
  end
  qs
end

function sort_window(a::AucBanditSelector)
  sort(as_vector(a.window), by = (index, reward) -> reward, rev = true)
end

# We should extend this to calculate all credits in one go instead of sorting 
# and looping for each index...
function area_under_curve_credit_assignment(a::AucBanditSelector, index)
  sorted_window = sort_window(a)
  wlength = length(sorted_window)
  y = q = 0.0
  rank = 1
  for (j, reward) in sorted_window
    # The AUC bandit skips the actual reward values and relies only on rank...
    reward = a.decay_factor^(rank-1) * (wlength - (rank - 1))
    if index == j
      y += reward
    else
      q += (y*reward)
    end
    rank += 1
  end
  q
end