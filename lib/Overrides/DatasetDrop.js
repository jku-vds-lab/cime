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
import { Grid } from "@mui/material";
import { CSVLoader, JSONLoader } from "projection-space-explorer";
import { DragAndDrop } from "projection-space-explorer";
import React from "react";
import { SDFLoader } from "./SDFLoader";
import { SDFModifierDialog } from "./SDFModifierDialog";
export var DatasetDrop = function (_a) {
    var onChange = _a.onChange, cancellablePromise = _a.cancellablePromise, abort_controller = _a.abort_controller;
    var _b = React.useState(null), entry = _b[0], setEntry = _b[1];
    var _c = React.useState(false), openSDFDialog = _c[0], setOpen = _c[1];
    function onModifierDialogClose(modifiers) {
        setOpen(false);
        if (modifiers !== null) {
            abort_controller = new AbortController();
            new SDFLoader().resolveContent(entry, onChange, cancellablePromise, modifiers, abort_controller);
        }
    }
    return (_jsxs(Grid, __assign({ container: true, item: true, alignItems: "stretch", justifyContent: "center", direction: "column", style: { padding: "16px" } }, { children: [_jsx(DragAndDrop, __assign({ accept: "*", handleDrop: function (files) {
                    if (files == null || files.length <= 0) {
                        return;
                    }
                    var file = files[0];
                    var fileName = file.name;
                    if (fileName.endsWith("sdf")) {
                        setEntry(file);
                        setOpen(true);
                    }
                    else {
                        var reader = new FileReader();
                        reader.onload = function (event) {
                            var _a;
                            var content = (_a = event === null || event === void 0 ? void 0 : event.target) === null || _a === void 0 ? void 0 : _a.result;
                            if (fileName.endsWith("json")) {
                                new JSONLoader().resolveContent(content, onChange);
                            }
                            else {
                                new CSVLoader().resolveContent(content, onChange);
                            }
                        };
                        reader.readAsText(file);
                    }
                } }, { children: _jsx("div", { style: { height: 200 } }, void 0) }), void 0), _jsx(SDFModifierDialog, { openSDFDialog: openSDFDialog, handleClose: onModifierDialogClose }, void 0)] }), void 0));
};
//# sourceMappingURL=DatasetDrop.js.map