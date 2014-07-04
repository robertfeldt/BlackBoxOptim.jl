#Julia = "julia"
Julia = "julia03"

desc "Run normal (fast) tests"
task :runtest do
  sh "#{Julia} -L src/BlackBoxOptim.jl test/runtests.jl"
end

desc "Run slow tests"
task :runslowtest do
  sh "#{Julia} -L src/BlackBoxOptim.jl test/runslowtests.jl"
end

desc "Run all tests"
task :runalltest => [:runtest, :runslowtest]

desc "Compare optimizers on standard, example problems"
task :compare_optimizers do
  sh "#{Julia} -L src/BlackBoxOptim.jl -L test/helper.jl test/test_compare_optimizers.jl"
end

def filter_latest_changed_files(filenames, numLatestChangedToInclude = 1)
  filenames.sort_by{ |f| File.mtime(f) }[-numLatestChangedToInclude, numLatestChangedToInclude]
end

desc "Run only the latest changed test file"
task :t do
  latest_changed_test_file = filter_latest_changed_files Dir["test/**/test*.jl"]
  sh "#{Julia} -L src/BlackBoxOptim.jl -L test/helper.jl #{latest_changed_test_file.first}"
end

desc "Clear build files etc"
task :clobber do
  Dir['**/*.jl.cov'].each do |f| 
    puts "Deleting #{f}"
    File.delete(f)
  end
end

task :at => :runalltest
task :st => :runslowtest

task :default => :runtest