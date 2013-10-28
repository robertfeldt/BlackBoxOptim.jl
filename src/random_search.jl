type RandomSearcher <: Optimizer
  search_space::Array{(Float64,Float64),1}
end

function ask(rs::RandomSearcher)
  
end
