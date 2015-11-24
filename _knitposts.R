#!/usr/bin/Rscript --vanilla

# Acknowledgement:
# Taken form ouzor.github.io


# compiles all .Rmd files in _R directory into .md files in _posts directory,
# if the input file is older than the output file.

# run ./knitpages.R to update all knitr files that need to be updated.


# path to Rmd files
# blog_path <- "blog/" # this works if you are in the "vladpetyuk.github.io"

# I'll use absolute for the same convenience of running from any directory.
# Perhaps I'll update this later to make sure it runs on multiple different comps.
blog_path <- "/Users/d3m629/Google Drive/vladpetyuk.github.io/blog/" # absolute!!



KnitPost <- function(input, outfile, base.url="/") {
  # this function is a modified version of an example here:
  # http://jfisher-usgs.github.com/r/2012/07/03/knitr-jekyll/
  require(knitr)
  opts_knit$set(base.dir = "/Users/d3m629/Google Drive/vladpetyuk.github.io")
  opts_knit$set(base.url = base.url)
  fig.path <- paste0("blog/figs/", sub(".Rmd$", "", basename(input)), "/")
  # fig.path <- paste0(blog_path, "figs/", sub(".Rmd$", "", basename(input)), "/")
  opts_chunk$set(fig.path = fig.path)
  opts_chunk$set(fig.cap = "figure")
  render_jekyll()
  knit(input, outfile, envir = parent.frame())
}


for (infile in list.files(paste0(blog_path,"_R/"), pattern="*.Rmd", full.names=TRUE)) {
  outfile = paste0(blog_path,"_posts/", sub(".Rmd$", ".md", basename(infile)))
  
  # knit only if the input file is the last one modified
  if (!file.exists(outfile) |
        file.info(infile)$mtime > file.info(outfile)$mtime) {
    KnitPost(infile, outfile)
  }
}

