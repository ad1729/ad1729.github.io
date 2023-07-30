
This project contains all the files for the blog.

## Instructions (for my future self)

- Create a folder under `posts/` with the folder name indicating the slug of the post
- Create an `ignored.qmd`[^1] file in the newly created folder
- Start writing
- When the post is ready for editing / proofreading, rename the file to `index.qmd`
- Run `quarto preview` in the terminal
- When satisfied with the result, send `Ctrl+C` to the terminal to cancel `quarto preview`
- Run styler on the source files (.R, .qmd)
- Run `quarto render`
- Commit changes and push to master / main branch
- Run `quarto publish gh-pages`
- Go do something else

[^1]: File is called `ignored.qmd` because of the `render: "!ignored.qmd"` parameter in the `_quarto.yml` file

## Extras

- Should the footnote be capitalized and end in a full stop? See [this](https://english.stackexchange.com/questions/242129/should-the-footnote-be-capitalized) thread.
