task :runtest do
  sh "julia -L src/GlobalOptim.jl test/runtests.jl"
end

task :runslowtest do
  sh "julia -L src/GlobalOptim.jl test/runslowtests.jl"
end

task :runalltest => [:runtest, :runslowtest]

task :t => :runtest
task :at => :runalltest
task :st => :runslowtest

task :default => :runtest