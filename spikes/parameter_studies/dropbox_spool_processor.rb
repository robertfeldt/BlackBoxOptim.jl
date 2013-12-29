# This is a simple spool processor for running julia optimization jobs using
# dropbox folders for managing the jobs.
require 'fileutils'

def spool_process(jobpath)
  basename = File.basename(jobpath)
  indir = File.dirname(jobname)
  basedir = File.dirname(indir)
  workdir = File.join(basedir, "work")
  outdir = File.join(basedir, "out")
  resultdir = File.join(basedir, "results")

  # Move the job to the work dir
  start_time = Time.now
  workname = start_time.strftime("%Y%m%d_%H%M%S") + "_#{machine_name}_" + basename
  workpath = File.join(workdir, workname)
  FileUtils.mv(jobpath, workpath)

  # Now run the job from the "results" dir
  Dir.cd(resultdir)
  system(workpath + " > " + workname + ".log")

  # Finished! We move the job to the out dir.
  end_time = Time.now
  outname = File.join(outdir, end_time.strftime("%Y%m%d_%H%M%S") + "_#{machine_name}_" + workname)
  outpath = File.join(outdir, outname)
  FileUtils.mv(workpath, outpath)

  puts "Job #{jobname} finished! elapsed = #{end_time-start_time}, started = #{start_time}"
end

SpoolDir = ARGV[1] || "/Users/feldt/Dropbox/job_processor"
WorkDir = File.join(SpoolDir, "work")

while(true)
  jobs = Dir.glob(WorkDir + "/*")
  if jobs.length < 1
    sleep(1.0 + 5.0 * rand())
  else
    begin
      job = jobs[rand(jobs.length)]
      spool_process(job)
    rescue => exception
      puts "Exception when processing job #{job}!"
      puts exception.backtrace
    end
  end
end
