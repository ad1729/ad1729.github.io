
This project contains all the files for the blog.

## Instructions (for my future self)

- Create a folder under `posts/` with the folder name indicating the slug of the post
- Create an `ignored.qmd` file in the newly created folder
  - File is called `ignored.qmd` because of the following parameters in the `_quarto.yml` file: 
  
  ```
  project:
  ...
  render:
    - "*.qmd"
    - "!ignored.qmd"
    - "!README.md"
    - "!LICENSE.md"
  ...
  ```
  
- Start writing
- When the post is ready for editing / proofreading, rename the file to `index.qmd`
- Run `quarto preview` in the terminal
- When satisfied with the result, send `Ctrl+C` to the terminal to cancel `quarto preview`
- Run `quarto render`
- Commit changes and push to master / main branch
- Run `quarto publish gh-pages`
- Go do something else

## Extras

- Should the footnote be capitalized and end in a full stop? See [this](https://english.stackexchange.com/questions/242129/should-the-footnote-be-capitalized) thread.
