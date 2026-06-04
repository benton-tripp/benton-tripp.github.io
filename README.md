# Benton Tripp's Portfolio & Blog

Check out my personal portfolio and blog here: [https://benton-tripp.github.io/](https://benton-tripp.github.io/)

## Local development

This Jekyll site uses Bundler with repo-local gems:

```powershell
bundle config set path vendor/bundle
bundle install
bundle exec jekyll serve
```

On Windows, make sure `ruby` and `bundle` resolve to a normal RubyInstaller/DevKit Ruby, not an embedded application Ruby. Check with:

```powershell
where ruby
where bundle
```

If these point to `C:\Program Files\ServiceNow\agent-client-collector\embedded\...`, move your RubyInstaller `bin` directory earlier on `PATH`, for example `C:\Ruby33-x64\bin`, then open a new terminal and rerun `bundle install`.
