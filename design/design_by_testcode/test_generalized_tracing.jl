using BlackBoxOptim
using DataFrames, CSV

# Proposal for more general tracing that can be used both to customize the
# current reporting, save to CSV files, and, longer-term, send information
# to GUI/plotting front-ends.
@testset "Generalized tracing interface" begin

@testset "Print trace to STDOUT (as per default) but also save to CSV file" begin
    # Set up a trace object that will print to STDOUT as well as save key data to CSV file.
    dtr = DefaultReportTracer(STDOUT, mode = :verbose)
    csvfilename = "mylog.csv"
    isfile(csvfilename) && rm(csvfilename; force = true) # Ensure file not exists
    csvtr = CsvFileTracer(csvfilename)
    tracer = SeqOfTracers(dtr, csvtr)

    # Now give the tracer we setup when setting up for optimization
    myfn(x) = sum(x)
    res = bboptimize(myfn; Tracer = tracer, SearchRange = (-1.0, 1.0), 
                NumDimensions = 2, MaxEvals = 100)
    @test isfile(csvfilename)
    df = CSV.read(csvfilename)
    @test size(df, 1) == 100 # As many rows as there are evaluations...
    @test df[end, :EvalNum] == 100
    @test df[end, :Fitness]     ~ 2.0 # Sum of 1.0 and 1.0
    @test df[end, :BestFitness] ~ 2.0 # Sum of 1.0 and 1.0
end

end