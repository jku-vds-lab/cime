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
import { Button, Grid, IconButton, List, ListItem, ListItemSecondaryAction, ListItemText, ListSubheader, } from "@mui/material";
import DeleteIcon from "@mui/icons-material/Delete";
import { DatasetType, useCancellablePromise } from "projection-space-explorer";
import React from "react";
import RefreshIcon from "@mui/icons-material/Refresh";
import { trackPromise } from "react-promise-tracker";
import { DEMO } from "../CIMEApp";
import { CIMEBackendFromEnv } from "../Backend/CIMEBackend";
import { LoadingIndicatorView } from "./DatasetTabPanel";
export var UploadedFiles = function (_a) {
    var onChange = _a.onChange, refresh = _a.refresh;
    var _b = React.useState([]), files = _b[0], setFiles = _b[1];
    var cancellablePromise = useCancellablePromise().cancellablePromise;
    React.useEffect(function () {
        updateFiles();
    }, [refresh]);
    var handleClick = function (entry) {
        onChange(entry);
    };
    var loading_area = "update_uploaded_files_list";
    function updateFiles() {
        trackPromise(cancellablePromise(CIMEBackendFromEnv.getFiles())
            .then(function (data) { return setFiles(data); })
            .catch(function (error) { return console.log(error); }), loading_area);
    }
    var handleDelete = function (file) {
        cancellablePromise(CIMEBackendFromEnv.deleteFile(file.id))
            .then(function () { return updateFiles(); })
            .catch(function (error) { return console.log(error); });
    };
    return (files && (_jsx("div", { children: _jsxs(Grid, __assign({ item: true, style: { overflowY: "auto", flex: "1 1 auto", maxHeight: "400px" } }, { children: [_jsxs(List, __assign({ subheader: _jsx("li", {}, void 0), style: { backgroundColor: "white" } }, { children: [!DEMO && (_jsxs(ListSubheader, { children: ["Uploaded Files", " ", _jsx(Button, __assign({ onClick: function () { return updateFiles(); } }, { children: _jsx(RefreshIcon, { style: { fontSize: "1.25rem" } }, void 0) }), void 0)] }, void 0)), DEMO && _jsx(ListSubheader, { children: "Select Dataset" }, void 0), files.map(function (file) { return (_jsxs(ListItem, __assign({ button: true, onClick: function () {
                                handleClick({
                                    display: file.name,
                                    path: file.id,
                                    type: DatasetType.Chem,
                                    uploaded: true, // indicates that file is already uploaded
                                });
                            } }, { children: [_jsx(ListItemText, { primary: file.name }, void 0), !DEMO && (_jsx(ListItemSecondaryAction, __assign({ onClick: function () {
                                        handleDelete(file);
                                    } }, { children: _jsx(IconButton, __assign({ edge: "end", "aria-label": "delete" }, { children: _jsx(DeleteIcon, {}, void 0) }), void 0) }), void 0))] }), file.id)); })] }), void 0), _jsx(LoadingIndicatorView, { area: loading_area }, void 0)] }), void 0) }, void 0))); //<Divider/>
};
//# sourceMappingURL=UploadedFiles.js.map