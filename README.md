# XYZeq: Spatially-resolved Single-cell RNA-sequencing

## Collaborators:

Seeing as how this repo is not public yet, I wanted to outline here how I think we should plan to make it ready for publication. Since we worked primarily in Python (and more specifically Jupyter Lab), I think this repo can just hold various notebooks each with the ability to reproduce the analysis for each figure in the paper. 

Since none of us have really been actively updating the notebooks here, I would suggest to clone this repo in a separate directory on your computer and work with it there. This avoids the risk of losing your work or overwriting everyone elses. Only add and edit files when you've tested with your own versions. I am **_not_** a Git expert, and welcome anyone who is to suggest better practices.

I realize we may be updating a lot of the analyses given the contamination we're trying to account for. I think before we do anything, we should make sure we can reproduce all of the figures in the current paper. Then, we can adjust the notebooks or add new ones to make new figures (if we feel it is necessary for resubmission).

In any case, I would like to propose the following organization. Please review:

### Data

#### Reproducing Figures

Raw and processed data will likely be deposited in some public genomics repository like GEO. Putting the raw data there is just a matter of uploading. However, processed data is less simple. We've done many analyses throughout our time working on the paper and we have many versions of the same data. To limit the number of files we put there, we should find out how to exactly reproduce every figure and only upload the data necessary for reproduction. 

#### Portability

In terms of portability, I think that once we upload the processed data, it should be a single tarball and therefore with a defined directory structure. When a person downloads that data, they simply have to change the path in /spatial/path.to.data.txt. What would be ideal is at the beginning of each notebook, the first cell defines a `path_prefix` variable that reads that path in and, since directory structure will be consistent, they can run all the notebooks without a problem, because all the file paths in each notebook should be relative to that.

### Notebooks

I've gone through all the notebooks and found the ones on my machine that are able to reproduce some of the figures. I've tried to clean those notebooks up, i.e. added markdown to explain what's happening, made figures a consistent size, got rid of empty cells, and got rid of superfluous cells with code that wasn't critical to reproduction of the figure. I Invite everyone to go through and revise the notebooks, add more details, etc. as long as the figures are reproduced consistently. 

I've placed those notebooks in /spatial/notebooks/figures for now. I haven't deleted anything just yet -- I put all the ones that were in there originally into their own directory (/spatial/notebooks/figures/original). The ones directly in /spatial/notebooks/figures/ are cleaned up. Let's keep everything there for now, and then when we're ready to make it public we can go throught the various other directories and delete unnecessary files and restructure.

### Environment

I'm not sure about you all, but I've been using a single conda environment for running all the analyses. I don't think I need to explain why environments are useful, but what might be worth including in this README is how to set conda environment as a Jupyter kernel, which is what I've been doing to run all my analyses. Of course, we should include a requirements file in this repo that is able to install all the packages (with specific version numbers) necessary to do all analyses.

### Final Testing

Right before we're ready to publish, I will open a new AWS machine, follow the instructions and install all the data and dependencies and run all the notebooks and make sure the figures look exactly like they do in the paper.

## Reproducing Figures:

All figure panels are in bold. Figure panels that have been reproduced successfully should be changed to italics. If you feel a certain panel doesn't really require a notebook to reproduce (e.g. it's just an image of the liver slice, or it was chiefly made in Illustrator), please also italicize it.

### 1
**_A (Illustrator)_**

**_B (Illustrator)_**

**C**

**D**

**E**

#### 2
**A**

**B**

**_C (image)_**

**D**

#### 3
**A**

**B**

**C**

**D**

**E**

**F**

#### 4
**A**

**B**

**C**

**D**
