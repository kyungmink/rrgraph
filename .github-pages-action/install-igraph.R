install.packages(
    'igraph', dependencies = c(
        'Depends', 'Imports', 'LinkingTo'
    ), INSTALL_opts = c(
        '--no-R', '--no-libs', '--no-data', '--no-help', '--no-demo',
        '--no-exec', '--no-inst', '--no-multiarch',
        '--without-keep.source', '--without-keep.parse.data', '--no-test-load',
        '--no-clean-on-error', '--no-staged-install'
    )
)
