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
import { Button, Dialog, DialogActions, DialogContent, DialogContentText, DialogTitle, TextField, } from "@mui/material";
import { connect } from "react-redux";
import React from "react";
import { setLineUpInput_dump } from "../State/LineUpInputDuck";
import { setDetailVisibility } from "projection-space-explorer";
var mapDispatchToProps = function (dispatch) { return ({
    setLineUp_dump: function (dump) { return dispatch(setLineUpInput_dump(dump)); },
    setLineUp_visibility: function (vis) { return dispatch(setDetailVisibility(vis)); },
}); };
var connector = connect(null, mapDispatchToProps);
// const LineUpContext = connector(function ({ lineUpInput, currentAggregation, setCurrentAggregation, setLineUpInput_visibility, onFilter, activeStory, hoverUpdate, hoverState }: Props)
export var LineUpDumpDialog = connector(function (_a) {
    var openDialog = _a.openDialog, setOpenDumpDialog = _a.setOpenDumpDialog, setLineUp_dump = _a.setLineUp_dump, setLineUp_visibility = _a.setLineUp_visibility;
    var _b = React.useState(""), dump = _b[0], setDump = _b[1];
    function handleChange(event) {
        setDump(event.target.value);
    }
    function handleClose() {
        setOpenDumpDialog(function () { return false; });
        setLineUp_visibility(true);
        setLineUp_dump(dump);
    }
    return (_jsxs(Dialog, __assign({ maxWidth: "lg", open: openDialog, onClose: handleClose }, { children: [_jsx(DialogTitle, { children: "Specify Modifiers" }, void 0), _jsxs(DialogContent, { children: [_jsx(DialogContentText, { children: "Insert Linup JSON dump" }, void 0), _jsx(TextField, { autoFocus: true, margin: "dense", id: "modifiers", label: "Modifiers", value: dump, onChange: handleChange, fullWidth: true }, void 0)] }, void 0), _jsxs(DialogActions, { children: [_jsx(Button, __assign({ onClick: handleClose }, { children: "Cancel" }), void 0), _jsx(Button, __assign({ onClick: handleClose }, { children: "Start" }), void 0)] }, void 0)] }), void 0));
});
//# sourceMappingURL=LineUpDumpDialog.js.map