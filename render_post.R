
library(readr)
library(magrittr)

render.post <- function(config.file="benton-tripp.github.io/_templates/config.yml",
                        template.file="benton-tripp.github.io/_templates/template.html",
                        input.rmd="benton-tripp.github.io/_Rmd/post.Rmd",
                        output.html="post.html",
                        output.dir = "benton-tripp.github.io/posts", 
                        knit.root.dir = "C:/users/bento/gis630",
                        title=NULL,
                        relative_path=NULL,
                        publish_time=format(Sys.time(), "%Y-%m-%dT%00:%00:00+00:00"),
                        publish_time_formatted = format(Sys.Date(), "%b %d, %Y"),
                        prev.relative_path=NULL,
                        prev.summary=NULL,
                        prev.title=NULL,
                        prev.description=NULL,
                        next.relative_path=NULL,
                        next.summary=NULL,
                        next.title=NULL,
                        next.description=NULL,
                        footer.p=NULL) {

    config <- yaml::read_yaml(config.file)
    
    if (!is.null(title)) config$title <- title
    if (!is.null(relative_path)) config$relative_path <- relative_path
    if (!is.null(publish_time)) config$publish_time <- publish_time
    if (!is.null(publish_time_formatted)) config$publish_time_formatted <- publish_time_formatted
    if (!is.null(prev.relative_path)) config$prev$relative_path <- prev.relative_path
    if (!is.null(prev.summary)) config$prev$summary <- prev.summary
    if (!is.null(prev.title)) config$prev$title <- prev.title
    if (!is.null(prev.description)) config$prev$description <- prev.description
    if (!is.null(next.relative_path)) config[["next"]]$relative_path <- next.relative_path
    if (!is.null(next.summary)) config[["next"]]$summary <- next.summary
    if (!is.null(next.title)) config[["next"]]$title <- next.title
    if (!is.null(next.description)) config[["next"]]$description <- next.description
    if (!is.null(footer.p)) config$footer$p <- footer.p
    
    template <- readr::read_file(template.file) %>%
      gsub("\\{\\{title\\}\\}", config$title, .) %>%
      gsub("\\{\\{relative_path\\}\\}", config$relative_path, .) %>%
      gsub("\\{\\{publish_time\\}\\}", config$publish_time, .) %>%
      gsub("\\{\\{publish_time_formatted\\}\\}", config$publish_time_formatted, .) %>%
      gsub("\\{\\{prev relative_path\\}\\}", config$prev$relative_path, .) %>%
      gsub("\\{\\{prev summary\\}\\}", config$prev$summary, .) %>%
      gsub("\\{\\{prev title\\}\\}", config$prev$title, .) %>%
      gsub("\\{\\{prev description\\}\\}", config$prev$description, .) %>%
      gsub("\\{\\{next relative_path\\}\\}", config[["next"]]$relative_path, .) %>%
      gsub("\\{\\{next summary\\}\\}", config[["next"]]$summary, .) %>%
      gsub("\\{\\{next title\\}\\}", config[["next"]]$title, .) %>%
      gsub("\\{\\{next description\\}\\}", config[["next"]]$description, .) %>%
      gsub("\\{\\{footer p\\}\\}", config$footer$p, .)
      
    
    cwd <- getwd()
    setwd(knit.root.dir)
    browser()
    rmarkdown::render(
      input = input.rmd, 
      output_format = "html_document", 
      output_file = output.html,
      output_dir = output.dir, 
      knit_root_dir = knit.root.dir,
      output_options = list(html_preview = F)
    )
    setwd(cwd)
    
    rendered.html.txt <- readr::read_file(file.path(output.dir, output.html))
    browser()
}

render.post(config.file="benton-tripp.github.io/_templates/config.yml",
            template.file="benton-tripp.github.io/_templates/template.html",
            input.rmd="C:/Users/bento/gis630/notebooks/eda.Rmd",
            output.html="2023-09-17-sdm-benchmark-study-part-2-exploratory-analysis.html",
            output.dir = "benton-tripp.github.io/posts", 
            knit.root.dir = "C:/users/bento/gis630",
            title="SDM Benchmark Study Part 2: Exploratory Analysis",
            relative_path="posts/2023-09-17-sdm-benchmark-study-part-2-exploratory-analysis.html",
            prev.relative_path="posts/sdm-benchmark-study-part-1-data-preparation.html",
            prev.title="SDM Benchmark Study Part 1: Data Preparation",
            prev.description="SDM Benchmark Study Part 1: Data Preparation")



