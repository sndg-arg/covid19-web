const path = require('path');

module.exports = {
    entry: './sndg_covid19/static/js/project.js',
    output: {
        filename: 'main.js',
        path: path.resolve(__dirname, 'staticfiles'),
    },
    node: {
        fs: 'empty'
    },
    module: {
        rules: [
            {
                test: /\.css$/i,
                use: ['style-loader', 'css-loader']
            }, {
                test: /\.(ttf|eot|svg|gif|png)(\?v=[0-9]\.[0-9]\.[0-9])?$/,
                use: [{
                    loader: 'file-loader'
                }]
            }
        ],
    },
    resolve: {
        alias: {
            'handlebars': 'handlebars/dist/handlebars.js'
            // Libraries msa and seqview are old and use require on handlebars that throws a warning
        }
    },
    mode:'development',
    devtool: 'source-map',
    stats: {
        warningsFilter: [
            'filter',
            /filter/,
            (warning) => warning.indexOf('msa/lib/') != -1
            // Library msa is old so it throws many unfixable warnings
        ]
    }
};
