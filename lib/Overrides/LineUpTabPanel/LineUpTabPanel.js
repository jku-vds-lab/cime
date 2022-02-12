var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { Typography } from "@mui/material";
import { Box, Button, FormControlLabel, Switch } from "@mui/material";
import { connect } from "react-redux";
import GetAppIcon from "@mui/icons-material/GetApp";
import { setLineUpInput_filter, } from "../../State/LineUpInputDuck";
import React from "react";
import { setDetailVisibility } from "projection-space-explorer";
var mapStateToProps = function (state) { return ({
    dataset: state.dataset,
    currentAggregation: state.currentAggregation,
    lineUpInput: state.lineUpInput,
}); };
var mapDispatchToProps = function (dispatch) { return ({
    setDetailVisibility: function (value) { return dispatch(setDetailVisibility(value)); },
    setLineUpInput_filter: function (value) { return dispatch(setLineUpInput_filter(value)); },
}); };
var connector = connect(mapStateToProps, mapDispatchToProps);
export var LineUpTabPanel = connector(function (_a) {
    var setDetailVisibility = _a.setDetailVisibility, setLineUpInput_filter = _a.setLineUpInput_filter, lineUpInput = _a.lineUpInput, dataset = _a.dataset, currentAggregation = _a.currentAggregation, splitRef = _a.splitRef;
    var handleChange = function (_, value) { };
    var onLoad = function (filter) {
        setDetailVisibility(true);
        setLineUpInput_filter(filter);
        var curr_sizes = splitRef.current.split.getSizes();
        if (curr_sizes[1] < 2) {
            splitRef.current.split.setSizes([curr_sizes[0], 70]);
        }
    };
    // https://stackoverflow.com/questions/31214677/download-a-reactjs-object-as-a-file
    var downloadImpl = function (data, name, mimetype) {
        var b = new Blob([data], { type: mimetype });
        var csvURL = window.URL.createObjectURL(b);
        var tempLink = document.createElement("a");
        tempLink.href = csvURL;
        tempLink.setAttribute("download", name);
        tempLink.click();
    };
    var exportCSV = function () {
        if (lineUpInput.lineup && lineUpInput.lineup.data) {
            // exports all data that is currently shown in the table -> filters and sorts are applied! also annotations are included
            lineUpInput
                .lineup.data.exportTable(lineUpInput.lineup.data.getRankings()[0], {
                separator: ",",
            })
                .then(function (x) { return downloadImpl(x, "lineup-export.csv", "application/csv"); });
        }
    };
    var _b = React.useState(false), cell_value_vis = _b[0], set_cell_value_vis = _b[1];
    React.useEffect(function () {
        var style = document.getElementById("cell_value_vis");
        if (!style) {
            style = document.createElement("style");
            style.setAttribute("id", "cell_value_vis");
            style.setAttribute("type", "text/css");
            var head = document.head || document.getElementsByTagName("head")[0];
            head.appendChild(style);
        }
        var css = cell_value_vis
            ? ".lu-hover-only { visibility: visible; }"
            : ".lu-hover-only { visibility: hidden; }";
        // @ts-ignore
        if (style.styleSheet) {
            // This is required for IE8 and below.
            // @ts-ignore
            style.styleSheet.cssText = css;
        }
        else {
            style.innerHTML = "";
            style.appendChild(document.createTextNode(css));
        }
    }, [cell_value_vis]);
    var toggleVis = function () {
        set_cell_value_vis(function () {
            return !cell_value_vis;
        });
    };
    return (_jsxs("div", __assign({ style: { display: "flex", flexDirection: "column", height: "100%" } }, { children: [_jsx(Box, __assign({ paddingLeft: 2, paddingTop: 1, paddingRight: 2 }, { children: _jsx(Typography, __assign({ variant: "subtitle2", gutterBottom: true }, { children: "LineUp Settings" }), void 0) }), void 0), _jsx(Box, __assign({ paddingLeft: 2, paddingTop: 1, paddingRight: 2 }, { children: _jsx(Button, __assign({ fullWidth: true, style: { marginRight: 2 }, variant: "outlined", onClick: function () { return onLoad({ reset: true }); } }, { children: "Load All" }), void 0) }), void 0), _jsx(Box, __assign({ paddingLeft: 2, paddingTop: 1, paddingRight: 2 }, { children: _jsx(Button, __assign({ fullWidth: true, variant: "outlined", onClick: function () {
                        return onLoad({ selection: currentAggregation.aggregation });
                    } }, { children: "Load Selection" }), void 0) }), void 0), _jsx(Box, __assign({ paddingLeft: 2, paddingTop: 1, paddingRight: 2 }, { children: _jsx(FormControlLabel, { control: _jsx(Switch, { color: "primary", value: cell_value_vis, onChange: function (event) {
                            toggleVis();
                        } }, void 0), label: "Show Cell Values" }, void 0) }), void 0), _jsx(Box, __assign({ paddingLeft: 2, paddingTop: 1, paddingRight: 2 }, { children: _jsxs(Button, __assign({ fullWidth: true, variant: "outlined", onClick: function () {
                        exportCSV();
                    } }, { children: [_jsx(GetAppIcon, {}, void 0), "\u00A0Export CSV"] }), void 0) }), void 0)] }), void 0));
});
function dispatch(arg0) {
    throw new Error("Function not implemented.");
}
//# sourceMappingURL=LineUpTabPanel.js.map