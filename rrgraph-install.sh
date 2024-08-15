#!/bin/bash
R CMD INSTALL --html rrgraph
sudo ln -f R/x86_64-pc-linux-gnu-library/3.3/rrgraph/html/* /var/www/html/man
sudo ln -f /var/www/html/man/00Index.html /var/www/html/man/index.html
ln -f rrgraph/R/* src
