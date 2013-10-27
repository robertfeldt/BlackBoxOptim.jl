task :runtest do
  sh "julia -L src/GlobalOptim.jl test/runtests.jl"
end

task :runslowtest do
  sh "julia -L src/GlobalOptim.jl test/runslowtests.jl"
end

task :runalltest => [:runtest, :runslowtest]

task :compare_optimizers do
  sh "julia -L src/GlobalOptim.jl test/compare_optimizers.jl"
end

def filter_latest_changed_files(filenames, numLatestChangedToInclude = 1)
  filenames.sort_by{ |f| File.mtime(f) }[-numLatestChangedToInclude, numLatestChangedToInclude]
end

desc "Run only the latest changed test file"
task :t do
  latest_changed_test_file = filter_latest_changed_files Dir["test/test*.jl"]
  sh "julia -L src/GlobalOptim.jl -L test/helper.jl #{latest_changed_test_file.first}"
end

task :at => :runalltest
task :st => :runslowtest

task :default => :runtest