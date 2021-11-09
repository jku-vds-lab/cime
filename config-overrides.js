/* global require, module, __dirname */
const { override, fixBabelImports, addWebpackAlias } = require('customize-cra')
const path = require('path')

module.exports = override(
  addWebpackAlias({
    'react': path.resolve('./node_modules/react'),
    'react-dom': path.resolve('./node_modules/react-dom'),
    'react-redux': path.resolve('./node_modules/react-redux'),
    '@mui/material': path.resolve('./node_modules/@mui/material'),
    '@emotion/styled': path.resolve('./node_modules/@emotion/styled'),
    '@emotion/react': path.resolve('./node_modules/@emotion/react'),
    'redux': path.resolve('./node_modules/redux'),
    '@mui/styles': path.resolve('./node_modules/@mui/styles'),
  })
)