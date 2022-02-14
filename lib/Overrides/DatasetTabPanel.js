import { Box, Button, Dialog, DialogActions, DialogContent, Typography, } from "@mui/material";
import { useCancellablePromise, RootActions } from "projection-space-explorer";
import React from "react";
import { usePromiseTracker } from "react-promise-tracker";
import { DatasetDrop } from "./DatasetDrop";
import { SDFLoader } from "./SDFLoader";
import { SDFModifierDialog } from "./SDFModifierDialog";
import { UploadedFiles } from "./UploadedFiles";
import Loader from "react-loader-spinner";
import { connect } from "react-redux";
export var LoadingIndicatorView = function (props) {
    var promiseInProgress = usePromiseTracker({ area: props.area }).promiseInProgress;
    return (promiseInProgress && (React.createElement("div", { style: {
            width: "100%",
            height: "100",
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
        } },
        React.createElement(Loader, { type: "ThreeDots", color: "#2BAD60", height: "100", width: "100" }))));
};
export var LoadingIndicatorDialog = function (props) {
    var promiseInProgress = usePromiseTracker({ area: props.area }).promiseInProgress;
    return (React.createElement(Dialog, { maxWidth: "lg", open: promiseInProgress },
        " ",
        React.createElement(DialogContent, null,
            React.createElement(LoadingIndicatorView, { area: props.area })),
        React.createElement(DialogActions, null,
            React.createElement(Button, { onClick: props.handleClose }, "Cancel"))));
};
var mapStateToProps = function (state) { return ({}); };
var mapDispatchToProps = function (dispatch) { return ({
    setDataset: function (value) { return dispatch(RootActions.loadDataset(value)); }
}); };
var connector = connect(mapStateToProps, mapDispatchToProps);
export var DatasetTabPanel = connector(function (_a) {
    var setDataset = _a.setDataset;
    var _b = React.useState(null), entry = _b[0], setEntry = _b[1];
    var _c = React.useState(false), openSDFDialog = _c[0], setOpen = _c[1];
    var _d = React.useState(0), refreshUploadedFiles = _d[0], setRefreshUploadedFiles = _d[1];
    var _e = useCancellablePromise(), cancellablePromise = _e.cancellablePromise, cancelPromises = _e.cancelPromises;
    var abort_controller = new AbortController();
    function onModifierDialogClose(modifiers) {
        //console.log(entry);
        setOpen(false);
        if (modifiers !== null) {
            abort_controller = new AbortController();
            new SDFLoader().resolvePath(entry, function (dataset) {
                //console.log("Dataset!!!", dataset);
                setDataset(dataset);
            }, cancellablePromise, modifiers, abort_controller);
        }
    }
    return (React.createElement("div", { style: { display: "flex", flexDirection: "column", height: "100%" } },
        React.createElement(Box, { paddingLeft: 2, paddingTop: 2 },
            React.createElement(Typography, { variant: "subtitle2", gutterBottom: true }, "Custom Datasets (Drag and Drop)")),
        React.createElement(DatasetDrop, { onChange: function (dataset) {
                //console.log("Dataaset", dataset);
                setDataset(dataset);
                setRefreshUploadedFiles(refreshUploadedFiles + 1);
            }, cancellablePromise: cancellablePromise, abort_controller: abort_controller }),
        React.createElement(Box, { paddingLeft: 2, paddingTop: 2 },
            React.createElement(Typography, { variant: "subtitle2", gutterBottom: true }, "Predefined Datasets")),
        React.createElement(UploadedFiles, { onChange: function (entry) {
                setEntry(entry);
                setOpen(true);
            }, refresh: refreshUploadedFiles }),
        React.createElement(LoadingIndicatorDialog, { handleClose: function () {
                cancelPromises();
            }, area: "global_loading_indicator" }),
        React.createElement(SDFModifierDialog, { openSDFDialog: openSDFDialog, handleClose: onModifierDialogClose })));
});
//# sourceMappingURL=DatasetTabPanel.js.map