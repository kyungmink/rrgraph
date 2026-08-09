install.packages(
    'package', repos = NULL, INSTALL_opts = c(
        '--html', '--no-R', '--no-libs', '--no-data', '--no-demo',
        '--no-exec', '--no-inst', '--no-multiarch',
        '--without-keep.source', '--without-keep.parse.data', '--no-test-load',
        '--no-clean-on-error', '--no-staged-install'
    )
)
file.rename(
    from = '/usr/local/lib/R/site-library/rrgraph/html/00Index.html',
    to   = '/usr/local/lib/R/site-library/rrgraph/html/index.html'
)
file.rename(
    from = '/usr/local/lib/R/site-library/rrgraph/help',
    to   = '/usr/local/lib/R/site-library/rrgraph/_help'
)
file.rename(
    from = '/usr/local/lib/R/site-library/rrgraph/html',
    to   = '/usr/local/lib/R/site-library/rrgraph/help'
)
file.copy(
    from = '/usr/local/lib/R/site-library/rrgraph/help',
    to   = '.', overwrite = TRUE, recursive = TRUE
)
