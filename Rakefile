Julia = "julia"
#Julia = "julia03"

Command = "#{Julia} --color=yes -L src/BlackBoxOptim.jl"

desc "Run normal (fast) tests"
task :runtest do
  sh "#{Command} test/runtests.jl"
end

desc "Run slow tests"
task :runslowtest do
  sh "#{Command} test/runslowtests.jl"
end

desc "Run all tests"
task :runalltest => [:runtest, :runslowtest]

desc "Compare optimizers on standard, example problems"
task :compare_optimizers do
  sh "#{Command} -L test/helper.jl test/test_compare_optimizers.jl"
end

def filter_latest_changed_files(filenames, numLatestChangedToInclude = 1)
  filenames.sort_by{ |f| File.mtime(f) }[-numLatestChangedToInclude, numLatestChangedToInclude]
end

desc "Run only the latest changed test file"
task :t do
  latest_changed_test_file = filter_latest_changed_files Dir["test/**/test*.jl"]
  sh "#{Command} -L test/helper.jl #{latest_changed_test_file.first}"
end

desc "Run and create code coverage information"
task :coverage do
  sh "#{Command} --code-coverage test/runtests.jl"
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