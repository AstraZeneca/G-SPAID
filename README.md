# G-SPAID
![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red) 


## Method 
The G-function computes the nearest neighbour distance for cells of type ‘j’ with respect to cells of type ‘i’ within a defined proximity distance (r). 
## Setting up the app
Download install R and R studio 
In R studio **click file** -> **New Project** -> **Version Control** -> **Git**  
Then paste link to github, add project directory name, and select the location where the app will downloaded and click **Create project**

- The app will be downloaded and console will read "The following package(s) will be installed: BiocManager [1.30.22] These packages will be installed into "/cell-cell-distribution/renv/library/R-4.1/x86_64-w64-mingw32"." type **y** and press **enter**
- type **renv::status()** and press enter to get a list of all the packages that will be installed
- type **rev::restore()** and press enter to install all the listed packages
- type **y** and enter to confirm installation 

After the installation is complete clich on the app.R file which will open the script 
click on the **Run App** button which will lauch the input window

## Data input 
1. **Choose CSV file** - HALO single cell data 
2. **Cell type:** - select column with cell types to perform cell-cell distribution. e.g. "F480 positive classification" and value greater than 0 used to identify postive cell type. 
3. **Image** - column to identify each sample
4. **group column** - if do not want to compare between 2 group, use the image column
5. **Xmin Xmax Ymin Ymax Column** - cells coordinate column
6. **distance** - value to calculate the G-area value. If used 30 the value would calculated from 0 - 30.
7. **no sim** - value indicating number of times to silmutate the distribution. To get alpha value of 0.05 (alpha = (2 ∗ 𝑛𝑟𝑎𝑛𝑘)/(1 + 𝑛𝑠𝑖𝑚)) use nsim = 199 and nrank = 5.
8. **rank** - m-th lowest and m-th highest. To get alpha value of 0.05 (alpha = (2 ∗ 𝑛𝑟𝑎𝑛𝑘)/(1 + 𝑛𝑠𝑖𝑚)) use nsim = 199 and nrank = 5. 

click **run** 
Once the concole prints "Done" close the input window

## Output
The app creates a folder starting with **Interaction** followed by date and time
The folder contains following files:
- **stat_table.csv** - file contains list of all interaction that are signifcantly different between the group with pvalue (wilcox test) and adjusted p-value (BH correction)  
- **heatmap.pdf** - pdf with all the cell-cell distribution analysed with y axis as 'i' cell type and x axis as 'j' cell type. Each i to j cell interaction is divided to show the scaled G-area value for both group with asterisk indicating padj signifance level.  
- Seperate folders for all significant interactions which contain:
  - **df_summary.csv** - file containing cell count of both i and j, G-area value, distance between the the observed curves is significant, average nearest neighbour distance for each sample in both group.
  - **g_*.pdf** - pdf file with g plot for each sample.
  - **h_*.pdf** - pdf contains density histogram for each sampl. 

