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
import { Box, Button, Dialog, DialogActions, DialogContent, Typography, } from "@mui/material";
import { useCancellablePromise } from "projection-space-explorer";
import React from "react";
import { usePromiseTracker } from "react-promise-tracker";
import { DatasetDrop } from "./DatasetDrop";
import { SDFLoader } from "./SDFLoader";
import { SDFModifierDialog } from "./SDFModifierDialog";
import { UploadedFiles } from "./UploadedFiles";
import Loader from "react-loader-spinner";
export var LoadingIndicatorView = function (props) {
    var promiseInProgress = usePromiseTracker({ area: props.area }).promiseInProgress;
    return (promiseInProgress && (_jsx("div", __assign({ style: {
            width: "100%",
            height: "100",
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
        } }, { children: _jsx(Loader, { type: "ThreeDots", color: "#2BAD60", height: "100", width: "100" }, void 0) }), void 0)));
};
export var LoadingIndicatorDialog = function (props) {
    var promiseInProgress = usePromiseTracker({ area: props.area }).promiseInProgress;
    return (_jsxs(Dialog, __assign({ maxWidth: "lg", open: promiseInProgress }, { children: [" ", _jsx(DialogContent, { children: _jsx(LoadingIndicatorView, { area: props.area }, void 0) }, void 0), _jsx(DialogActions, { children: _jsx(Button, __assign({ onClick: props.handleClose }, { children: "Cancel" }), void 0) }, void 0)] }), void 0));
};
export function DatasetTabPanel(_a) {
    var onDataSelected = _a.onDataSelected;
    var _b = React.useState(null), entry = _b[0], setEntry = _b[1];
    var _c = React.useState(false), openSDFDialog = _c[0], setOpen = _c[1];
    var _d = React.useState(0), refreshUploadedFiles = _d[0], setRefreshUploadedFiles = _d[1];
    var _e = useCancellablePromise(), cancellablePromise = _e.cancellablePromise, cancelPromises = _e.cancelPromises;
    var abort_controller = new AbortController();
    function onModifierDialogClose(modifiers) {
        console.log(entry);
        setOpen(false);
        if (modifiers !== null) {
            abort_controller = new AbortController();
            new SDFLoader().resolvePath(entry, function (dataset) {
                console.log("Dataset!!!", dataset);
                onDataSelected(dataset);
            }, cancellablePromise, modifiers, abort_controller);
        }
    }
    return (_jsxs("div", __assign({ style: { display: "flex", flexDirection: "column", height: "100%" } }, { children: [_jsx(Box, __assign({ paddingLeft: 2, paddingTop: 2 }, { children: _jsx(Typography, __assign({ variant: "subtitle2", gutterBottom: true }, { children: "Custom Datasets (Drag and Drop)" }), void 0) }), void 0), _jsx(DatasetDrop, { onChange: function (dataset) {
                    console.log("Dataaset", dataset);
                    onDataSelected(dataset);
                    setRefreshUploadedFiles(refreshUploadedFiles + 1);
                }, cancellablePromise: cancellablePromise, abort_controller: abort_controller }, void 0), _jsx(Box, __assign({ paddingLeft: 2, paddingTop: 2 }, { children: _jsx(Typography, __assign({ variant: "subtitle2", gutterBottom: true }, { children: "Predefined Datasets" }), void 0) }), void 0), _jsx(UploadedFiles, { onChange: function (entry) {
                    console.log(entry);
                    setEntry(entry);
                    setOpen(true);
                }, refresh: refreshUploadedFiles }, void 0), _jsx(LoadingIndicatorDialog, { handleClose: function () {
                    cancelPromises();
                }, area: "global_loading_indicator" }, void 0), _jsx(SDFModifierDialog, { openSDFDialog: openSDFDialog, handleClose: onModifierDialogClose }, void 0)] }), void 0));
}
//# sourceMappingURL=DatasetTabPanel.js.map