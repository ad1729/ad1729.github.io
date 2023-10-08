This project contains all the files for the blog.

## Instructions (for my future self)

-   Create a folder under `posts/` with the folder name indicating the slug of the post

-   Create an `ignored.qmd`[^readme-1] file in the newly created folder

    -   If using Python + Jupyter notebook from within VS Code, then create `ignored.ipynb` instead of `ignored.qmd`

-   Start writing

    -   If the post contains `Python` code, then the `venv` for the blog should be used. This is done automatically when using Quarto with RStudio (on my linux machine) and when using VS Code. See the section on managing python environments below.

-   For R - When the post is ready for editing / proofreading, rename the file to `index.qmd`

-   Run `quarto preview` in the terminal, or click the Render button and see the post in the Viewer pane within RStudio. Same works for VS Code.

-   When satisfied with the result, send `Ctrl+C` to the terminal to cancel `quarto preview`

-   For converting the final ipynb file, this command can be used from the terminal

    ```{quarto convert posts/<post_dir>/ignored.ipynb --output posts/<post_dir>/index.qmd}
    ```

-   Optionally run styler on the source files (.R, .qmd). Useful to stage the source file to see what styler has changed and whether it should be accepted or rejected.

-   Run `quarto render` if `quarto preview` was used to review the file

-   Commit changes and push to master / main branch

-   Run `quarto publish gh-pages`

-   Go do something else

[^readme-1]: File is called `ignored.qmd` because of the `render: "!ignored.qmd"` parameter in the `_quarto.yml` file

## Managing R environments

For this, the renv documentation can be referred to. Important to periodically run

```{r}
renv::status()
```

and

```{r}
renv::snapshot()
```

## Managing Python environments

Python packages are managed using pip. The local python environment for the blog can be activated by running the following from the root directory of the blog via bash

``` bash
source blogpyenv/bin/activate
```

See [this](https://docs.python.org/3/tutorial/venv.html) link for managing venv.

If new packages are installed or upgraded, then the `requirements.txt` file should be updated by running

``` bash
python -m pip freeze > requirements.txt
```

I should probably use poetry -- which I would for complex projects -- but will switch in the future if I end up having enough Python posts.

The `reticulate` package in R is key when writing and running Python code from within RStudio.

### Python with VS Code

See these links:

-   <https://quarto.org/docs/computations/python.html>

-   <https://quarto.org/docs/tools/jupyter-lab.html>

-   <https://quarto.org/docs/tools/vscode.html>

-   <https://quarto.org/docs/tools/vscode-notebook.html>

-   <https://quarto.org/docs/visual-editor/vscode/>

## Extras

-   Should the footnote be capitalized and end in a full stop? See [this](https://english.stackexchange.com/questions/242129/should-the-footnote-be-capitalized) thread.
