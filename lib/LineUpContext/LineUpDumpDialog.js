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
    return (React.createElement(Dialog, { maxWidth: "lg", open: openDialog, onClose: handleClose },
        React.createElement(DialogTitle, null, "Specify Modifiers"),
        React.createElement(DialogContent, null,
            React.createElement(DialogContentText, null, "Insert Linup JSON dump"),
            React.createElement(TextField, { autoFocus: true, margin: "dense", id: "modifiers", label: "Modifiers", value: dump, onChange: handleChange, fullWidth: true })),
        React.createElement(DialogActions, null,
            React.createElement(Button, { onClick: handleClose }, "Cancel"),
            React.createElement(Button, { onClick: handleClose }, "Start"))));
});
//# sourceMappingURL=LineUpDumpDialog.js.map