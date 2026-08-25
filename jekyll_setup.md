# Setting up Jekyll for Serving the Site

1. **Install Ruby:**  
   Jekyll is written in Ruby, so you need a recent version of Ruby installed on your system. You can download Ruby from [ruby-lang.org](https://www.ruby-lang.org/) or use a version manager like [rbenv](https://github.com/rbenv/rbenv) or [RVM](https://rvm.io/).

   You might need to do something like this to add to path:

   ```cmd
   set "PATH=C:\Ruby40-x64\bin;%PATH%"
   ```

2. **Install Bundler and Jekyll:**  
   Once Ruby is installed, open your terminal (or command prompt) and run:  
   
   ```bash
   gem install bundler jekyll
   ```

3. **Create a New Jekyll Site:**  
   Navigate to the directory where you want your site and run:  

   ```bash
   jekyll new portfolio-site
   cd portfolio-site
   ```

   This command creates a new directory named `portfolios-site` with a default structure (including a `_layouts` folder, `_posts` folder, etc.).

4. **Customize Your Layout:**  
   - Open the `_layouts/default.html` file to edit your overall template (including the `<head>`, navbar, footer, etc.).
   - Use Liquid templating syntax (e.g., `{{ content }}`) to designate where page-specific content will go.

5. **Create Pages:**  
   Create individual pages (for example, `index.md`, `about.md`, etc.) in the root directory. At the top of each file, add a front matter block:
   ```yaml
   ---
   layout: default
   title: Home
   ---
   ```
   Then add your page content below.

6. **Build and Serve Your Site:**  
   From the project directory, run:
   ```bash
   bundle exec jekyll serve
   ```
   This will build your site and serve it locally (typically at `http://localhost:4000`).

7. **Deployment:**  
   Once you’re happy with your site, you can deploy the static files (found in the `_site` directory after a build) to any static hosting provider (e.g., GitHub Pages, Netlify, Vercel).

### Summary

- **Installation:** Ruby → Bundler & Jekyll via `gem install`
- **Site Creation:** `jekyll new portfolio-site`
- **Customization:** Edit `_layouts/default.html` and add pages with front matter.
- **Local Development:** `bundle exec jekyll serve`
- **Deployment:** Upload the static files from the `_site` folder.

This setup produces full, SEO-friendly static HTML pages while allowing you to manage your layout with a templating language. If you prefer another static site generator (like Hugo, Eleventy, etc.), the setup process is similar—install the tool, create a new project, and customize templates.



