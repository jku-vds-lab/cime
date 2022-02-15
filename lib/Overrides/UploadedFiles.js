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
            .then(function (data) {
            setFiles(data !== null && data !== void 0 ? data : []);
        })
            .catch(function (error) { return console.log(error); }), loading_area);
    }
    var handleDelete = function (file) {
        cancellablePromise(CIMEBackendFromEnv.deleteFile(file.id))
            .then(function () { return updateFiles(); })
            .catch(function (error) { return console.log(error); });
    };
    return (files && (React.createElement("div", null,
        React.createElement(Grid, { item: true, style: { overflowY: "auto", flex: "1 1 auto", maxHeight: "400px" } },
            React.createElement(List, { subheader: React.createElement("li", null), style: { backgroundColor: "white" } },
                !DEMO && (React.createElement(ListSubheader, null,
                    "Uploaded Files",
                    " ",
                    React.createElement(Button, { onClick: function () { return updateFiles(); } },
                        React.createElement(RefreshIcon, { style: { fontSize: "1.25rem" } })))),
                DEMO && React.createElement(ListSubheader, null, "Select Dataset"),
                files.map(function (file) { return (React.createElement(ListItem, { key: file.id, button: true, onClick: function () {
                        handleClick({
                            display: file.name,
                            path: file.id,
                            type: DatasetType.Chem,
                            uploaded: true,
                        });
                    } },
                    React.createElement(ListItemText, { primary: file.name }),
                    !DEMO && (React.createElement(ListItemSecondaryAction, { onClick: function () {
                            handleDelete(file);
                        } },
                        React.createElement(IconButton, { edge: "end", "aria-label": "delete" },
                            React.createElement(DeleteIcon, null)))))); })),
            React.createElement(LoadingIndicatorView, { area: loading_area })))));
};
//# sourceMappingURL=UploadedFiles.js.map