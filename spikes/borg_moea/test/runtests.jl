using Autotest

# We need this import since there are conflicting describe methods in DataFrames and StatsBase...
import Autotest: describe

if length(ARGS) > 0 && ARGS[1] == "continuous"
  Autotest.autorun("PerfTest", "src", "test")
else
  Autotest.run("PerfTest", "src", "test")
end