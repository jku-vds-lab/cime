import { Button, Dialog, DialogActions, DialogContent, DialogContentText, DialogTitle, TextField, } from "@mui/material";
import React from "react";
export function SDFModifierDialog(_a) {
    var openSDFDialog = _a.openSDFDialog, handleClose = _a.handleClose;
    var _b = React.useState(""), modifiers = _b[0], setModifiers = _b[1];
    function handleModifierChange(event) {
        setModifiers(event.target.value);
    }
    return (React.createElement(Dialog, { maxWidth: "lg", open: openSDFDialog, onClose: function () { return handleClose(null); } },
        React.createElement(DialogTitle, null, "Specify Modifiers"),
        React.createElement(DialogContent, null,
            React.createElement(DialogContentText, null,
                "Manually specify modifiers separated by semicolons e.g. \"pred;fp;latent\". ",
                React.createElement("br", null),
                "You can also leave this field empty, if the modifiers are included by default. ",
                React.createElement("br", null),
                "The following modifiers are included by default: \"pred\", \"predicted\", \"measured\", \"fingerprint\", \"rep\"."),
            React.createElement(TextField, { autoFocus: true, margin: "dense", id: "modifiers", label: "Modifiers", value: modifiers, onChange: handleModifierChange, fullWidth: true })),
        React.createElement(DialogActions, null,
            React.createElement(Button, { onClick: function () { return handleClose(null); } }, "Cancel"),
            React.createElement(Button, { onClick: function () { return handleClose(modifiers); } }, "Start"))));
}
//# sourceMappingURL=SDFModifierDialog.js.map