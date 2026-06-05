# Benton Tripp's Portfolio & Blog

Check out my personal portfolio and blog here: [https://benton-tripp.github.io/](https://benton-tripp.github.io/)

## Local development

This Jekyll site uses Bundler with repo-local gems:

```cmd
bundle config set path vendor/bundle
bundle install
bundle exec jekyll serve
```

On Windows, make sure `ruby` and `bundle` resolve to a normal RubyInstaller/DevKit Ruby, not an embedded application Ruby. Check with:

```cmd
where ruby
where bundle
```

If these only point to `C:\Program Files\APPLICATION\...`, install RubyInstaller with DevKit from <https://rubyinstaller.org/downloads/> first. After installing, open a new terminal and confirm a RubyInstaller path such as `C:\Ruby40-x64\bin` appears before the ServiceNow path:

```cmd
where ruby
ruby -v
```

If needed, put RubyInstaller first for the current terminal session.

```cmd
set "PATH=C:\Ruby40-x64\bin;%PATH%"
```

Then install the locked Bundler version and start Jekyll:

```cmd
gem install bundler -v 2.6.3
bundle install
bundle exec jekyll serve
```

Ruby 4 may print Bundler warnings about already initialized `Gem::Platform` constants. These are harmless if Jekyll continues and reports a successful build or serve.
