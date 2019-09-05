
# setup info here:
# https://askubuntu.com/questions/218708/installing-latest-version-of-r-base
# and updating apt repository + gpg key from here:
# https://cran.r-project.org/bin/linux/ubuntu/README.html
sudo apt update
add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial-cran35/"
gpg --keyserver keyserver.ubuntu.com --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -

sudo apt-key update
sudo apt-get update

sudo apt install r-base r-base-dev gdebi-core libxml2-dev libssl-dev gcc g++ libcurl4-openssl-dev

wget https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.1335-amd64.deb

sudo gdebi rstudio-server-1.2.1335-amd64.deb
rm rstudio-server-1.2.1335-amd64.deb

Rscript --vanilla install_r_packages.R
