# us-mortality
A project exploring US mortality by state using structural time series models

Email for working paper: paige_park@berkeley.edu

## Reproducing the analysis

This project uses [Pixi](https://pixi.sh) to manage R and system dependencies, and [renv](https://rstudio.github.io/renv/) to manage R packages.

### Setup

1. **Install Pixi** (if you don't have it):

   ```sh
   curl -fsSL https://pixi.sh/install.sh | sh
   ```

2. **Clone this repo and enter the directory:**

   ```sh
   git clone <your-repo-url>
   cd us-mortality
   ```

3. **Install the environment and R packages:**

   ```sh
   pixi install
   pixi run setup
   ```

That's it! You're ready to run the analysis.
