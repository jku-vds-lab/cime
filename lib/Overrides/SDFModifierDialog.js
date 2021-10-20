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
import React from "react";
export function SDFModifierDialog(_a) {
    var openSDFDialog = _a.openSDFDialog, handleClose = _a.handleClose;
    var _b = React.useState(""), modifiers = _b[0], setModifiers = _b[1];
    function handleModifierChange(event) {
        setModifiers(event.target.value);
    }
    return (_jsxs(Dialog, __assign({ maxWidth: "lg", open: openSDFDialog, onClose: function () { return handleClose(null); } }, { children: [_jsx(DialogTitle, { children: "Specify Modifiers" }, void 0), _jsxs(DialogContent, { children: [_jsxs(DialogContentText, { children: ["Manually specify modifiers separated by semicolons e.g. \"pred;fp;latent\". ", _jsx("br", {}, void 0), "You can also leave this field empty, if the modifiers are included by default. ", _jsx("br", {}, void 0), "The following modifiers are included by default: \"pred\", \"predicted\", \"measured\", \"fingerprint\", \"rep\"."] }, void 0), _jsx(TextField, { autoFocus: true, margin: "dense", id: "modifiers", label: "Modifiers", value: modifiers, onChange: handleModifierChange, fullWidth: true }, void 0)] }, void 0), _jsxs(DialogActions, { children: [_jsx(Button, __assign({ onClick: function () { return handleClose(null); } }, { children: "Cancel" }), void 0), _jsx(Button, __assign({ onClick: function () { return handleClose(modifiers); } }, { children: "Start" }), void 0)] }, void 0)] }), void 0));
}
//# sourceMappingURL=SDFModifierDialog.js.map