import React = require("react");
import { Typography, Slider } from "@material-ui/core";
import './SizeSlider.scss'
import { connect } from 'react-redux'
import { setGlobalPointSize } from "../../../Ducks/GlobalPointSizeDuck";

const SizeSliderFull = ({ sizeScale, setRange }) => {
    const marks = [
        {
            value: 0,
            label: '0',
        },
        {
            value: 1,
            label: `1`,
        },
        {
            value: 2,
            label: `2`,
        },
        {
            value: 3,
            label: `3`,
        },
        {
            value: 4,
            label: `4`,
        },
        {
            value: 5,
            label: `5`,
        },
    ];

    return <div className="SizeSliderParent">
        <Typography id="range-slider" gutterBottom>
            Size Scale
      </Typography>
        <Slider
            min={0}
            max={5}
            value={sizeScale}
            onChange={(_, newValue) => setRange(newValue)}
            step={0.25}
            marks={marks}
            valueLabelDisplay="auto"
        ></Slider>
    </div>
}

/**
 * 
 * @param state                 sizeScale={this.state.vectorBySize.values.range}
                onChange={(e, newVal) => {
                  if (arraysEqual(newVal, this.state.vectorBySize.values.range)) {
                    return;
                  }

                  this.state.vectorBySize.values.range = newVal

                  this.setState({
                    vectorBySize: this.state.vectorBySize
                  })

                  if (this.state.vectorBySize != null) {
                    this.threeRef.current.particles.sizeCat(this.state.vectorBySize)
                    this.threeRef.current.particles.updateSize()
                  }
                }}
 */

const mapStateToProps = state => ({
    sizeScale: state.globalPointSize
})

const mapDispatchToProps = dispatch => ({
    setRange: value => dispatch(setGlobalPointSize(value))
})

export const SizeSlider = connect(mapStateToProps, mapDispatchToProps)(SizeSliderFull)