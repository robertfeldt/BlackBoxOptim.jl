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

@testset "Trace also to a front-end/GUI that shows fitness progress etc" begin
    # Not sure this is the right design but since we need to update a front-end
    # basically on every trace event we could use a special tracer that runs both
    # a http server to serve the GUI and a websocket that communicates with that
    # http server to update the GUI. We could use any GUI, really, but a simple
    # one might be to put up a Vega-Lite graph that plots fitness over time
    # on logarithmic scales.

    dtr = DefaultReportTracer(STDOUT, mode = :verbose)
    webguitr = VegaLiteOptimizationFrontEnd(; port = 8082) # Open localhost:8082 to access GUI once the optimization has started
    tracer = SeqOfTracers(dtr, webguitr)
    myfn(x) = sum(x)
    res = bboptimize(myfn; Tracer = tracer, SearchRange = (-1.0, 1.0), 
                NumDimensions = 2, MaxEvals = 100)
    # We need to get the page and make some sanity checks on it
    try
        r = HTTP.request("GET", "http://localhost:8082"; readtimeout = 10)
        # Add some sanity checks of the gui html here...
    catch err
        @test false # Indicate that exception was thrown => failure
    end
end

end