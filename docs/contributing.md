# Editing This Site

This page gives documentation about how to make an edit to this site.

## 1) Install `git`
If it's not already installed on your system, the first step is to install the `git` utility. This is what underlies github.com and 

!!! Note
    Make sure to make a github account so that you can push 
    changes to public repositories

## 2) Clone this Repository

On your system, open up a folder and change the directory (`cd` command) to wherever you would like this repository to be saved. Then clone this repository and `cd` into the directory!

```shell
cd /home/smazurchuk/Desktop
git clone https://github.com/mcw-mstp/mcw-mstp.github.io.git
cd mcw-mstp.github.io
```

## 3) Create a New Working Branch

We now have to create a new banch. You can pick whatever name you like for this, and create the branch with the command

```shell
git checkout -b branch-name-here
```

The `-b` option indicates that you are creating a new branch.

## 4) Install mkdocs (Optional)

!!! info inline end
    Installing this is only necessary to preview changes

In order to generate a live preview of the edits you make to the site, you will need to download a python package called mkdocs. If python or a way of managing python virtual environments is not already installed on your system, there is better documentation **elsewhere** for doing those steps. 

If you have python installed, then simply type:

```shell
pip install mkdocs
```

To see if all the required packages installed, try to serve the webpage
```shell
mkdocs serve
```

If there is an error, there might be some extenstions missing. In this case, try the following command

```shell
pip install mkdocs-material-extensions
```

If not, shoot me an email

## 5) Make an edit!

We suggest that the first edit is adding your name to the [contributers](/people.md) page! To do this, open the file `people.md` file in the docs folder in your favorite text editor. Then add your name using | as a dividers between table values.

```pre
| Name | Email | Background |
|---|---|---|
| Stephen | smazurchuk@mcw.edu | I work mainly with fMRI datasets (imagingn data) and have a strong background in python and Matlab. <br>Also somewhat knowledgeable in some machine learning topics |
|Dohn Joe  | dj@mcw.edu | A made-up personality |
```

!!! Note
    If you would like for the text in the background section to be broken up over multiple lines,
    make sure to include the `<br>` text where you would like a line break

## 6) Push Changes to Github

This section will perhaps have the *strangest* instructions for those new to `git`. Now that you have edited and saved your changes to `people.md`. The steps are:

1. `git add docs/people.md`
2. `git commit -m "your message here in qoutes"`
3. `git push origin branch-name-here`

## 7) Make Pull Request for Official Site

Finally, this last step submits your changes for review before it is *merged* with the main branch. To this, you can go to github.com and when you go to this repository, you should see the ability to click a green button that says "Compare & Pull Request". Follow the steps, and write a description of what your pull request changes on the site. This info will be reviewed and then merged with the site! That's it! 

# Extra Info

Mkdocs is a very convienent tool for building documentation. For working with this site, the only commands you should need are:

* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build --force` - Build the documentation site. (This pushes the built site to the gh-pages branch). The force option is needed since we hace multiple pages under some of the main headers
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

A comprehensive list of features that can be used on this site can be found at: 

<https://squidfunk.github.io/mkdocs-material/reference/abbreviations/>

*[text editor]: Which of course is vim