var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
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
var __spreadArray = (this && this.__spreadArray) || function (to, from, pack) {
    if (pack || arguments.length === 2) for (var i = 0, l = from.length, ar; i < l; i++) {
        if (ar || !(i in from)) {
            if (!ar) ar = Array.prototype.slice.call(from, 0, i);
            ar[i] = from[i];
        }
    }
    return to.concat(ar || Array.prototype.slice.call(from));
};
import * as React from "react";
import "./chem.scss";
import { Autocomplete, Box, Button, Checkbox, createFilterOptions, FormControlLabel, Grid, IconButton, Paper, Popover, Slider, Switch, TextField, Tooltip, Typography, } from "@mui/material";
import { trackPromise } from "react-promise-tracker";
import { connect } from "react-redux";
import RefreshIcon from "@mui/icons-material/Refresh";
import SettingsIcon from "@mui/icons-material/Settings";
import InfoIcon from "@mui/icons-material/Info";
import FilterListIcon from "@mui/icons-material/FilterList";
import { isFunction } from "lodash";
import { useCancellablePromise } from "projection-space-explorer";
import AddCircleOutlineIcon from "@mui/icons-material/AddCircleOutline";
import DeleteIcon from "@mui/icons-material/Delete";
import { selectVectors, setHoverState, WindowMode, } from "projection-space-explorer";
import { CIMEBackendFromEnv } from "../../Backend/CIMEBackend";
import { LoadingIndicatorView } from "../../Overrides/DatasetTabPanel";
import { setRDKit_contourLines, setRDKit_doAlignment, setRDKit_refresh, setRDKit_scale, setRDKit_showMCS, setRDKit_sigma, setRDKit_width, } from "../../State/RDKitSettingsDuck";
/**
 * Chem Legend, implemented
 */
var mapStateToProps_Chem = function (state) {
    var _a;
    return ({
        dataset: state.dataset,
        hoverSettings: state.hoverSettings,
        rdkitSettings: state.rdkitSettings,
        columns: (_a = state.dataset) === null || _a === void 0 ? void 0 : _a.columns,
    });
};
var mapDispatchToProps_Chem = function (dispatch) { return ({
    setCurrentAggregation: function (samples) {
        return dispatch(selectVectors(samples, false));
    },
    setHoverstate: function (state, updater) {
        return dispatch(setHoverState(state, updater));
    },
}); };
var connector_Chem = connect(mapStateToProps_Chem, mapDispatchToProps_Chem);
export var ChemLegendParent = connector_Chem(function (props) {
    var _a = useCancellablePromise(), cancellablePromise = _a.cancellablePromise, cancelPromises = _a.cancelPromises;
    if (props.mcs_only) {
        var _b = React.useState(React.createElement("div", null, "loading...")), mcsComp = _b[0], setMcsComp_1 = _b[1];
        var smiles_col_1 = get_smiles_col(props.columns);
        React.useEffect(function () {
            cancelPromises();
            if (smiles_col_1 in props.columns) {
                var controller = new AbortController();
                var my_fetch = void 0;
                if (props.diff && props.selection_ref) {
                    var smilesA = props.selection.map(function (row) { return row[smiles_col_1]; });
                    var smilesB = props.selection_ref.map(function (row) { return row[smiles_col_1]; });
                    my_fetch = CIMEBackendFromEnv.getDifferenceHighlight(smilesA, smilesB, controller);
                }
                else {
                    var formData_1 = new FormData();
                    props.selection.every(function (row) {
                        formData_1.append("smiles_list", row[smiles_col_1]);
                        return true;
                    });
                    my_fetch = CIMEBackendFromEnv.getMCSFromSmilesList(formData_1, controller);
                }
                trackPromise(cancellablePromise(my_fetch.then(function (x) {
                    if (x.length > 100) {
                        // check if it is actually long enogh to be an img
                        setMcsComp_1(function () { return (React.createElement("div", { style: {
                                width: 200,
                                height: 200,
                                backgroundSize: "contain",
                                backgroundPosition: "center",
                                backgroundRepeat: "no-repeat",
                                backgroundImage: "url('data:image/jpg;base64,".concat(x, "')"),
                            } })); });
                    }
                    else {
                        setMcsComp_1(function () { return React.createElement("div", null, x); });
                    }
                }), controller));
            }
        }, [props.selection, props.selection_ref, props.mcs_only]);
        if (smiles_col_1 in props.columns) {
            return React.createElement("div", null, mcsComp);
        }
        return React.createElement("div", null, "No SMILES column found.");
    }
    var _c = React.useState(false), settingsOpen = _c[0], setSettingsOpen = _c[1];
    var _d = React.useState(["Common Substructure"]), repList = _d[0], setRepList = _d[1];
    var _e = React.useState([0]), chemComponents = _e[0], setChemComponents = _e[1];
    var _f = React.useState(["Common Substructure"]), chemComponentsCurrentRep = _f[0], setChemComponentsCurrentRep = _f[1];
    var loadRepList = function (refresh) {
        if (refresh === void 0) { refresh = false; }
        if (refresh || repList.length <= 1) {
            var loading_area_1 = "global_loading_indicator";
            var controller = new AbortController();
            trackPromise(cancellablePromise(CIMEBackendFromEnv.getRepresentationList(refresh, props.dataset.info.path, controller).then(function (x) {
                if (x["rep_list"].length > 0) {
                    var rep_list = __spreadArray([], x["rep_list"], true);
                    rep_list.splice(0, 0, "Common Substructure");
                    setRepList(rep_list);
                }
            }), controller), loading_area_1);
        }
    };
    var addComp = function () {
        var comps = __spreadArray([], chemComponents, true);
        comps.push(Math.max.apply(Math, comps) + 1);
        setChemComponents(comps);
        var compsCR = __spreadArray([], chemComponentsCurrentRep, true);
        compsCR.push("Common Substructure");
        setChemComponentsCurrentRep(compsCR);
    };
    React.useEffect(function () {
        cancelPromises();
        if (props.aggregate) {
            loadRepList();
        }
    }, []);
    var removeComponent = function (id) {
        var comps = __spreadArray([], chemComponents, true);
        var compsCR = __spreadArray([], chemComponentsCurrentRep, true);
        var index = comps.indexOf(id);
        if (index > -1) {
            comps.splice(index, 1);
            compsCR.splice(index, 1);
        }
        setChemComponents(comps);
        setChemComponentsCurrentRep(compsCR);
    };
    var setCurrentRep = function (value, id) {
        if (repList.includes(value)) {
            var compsCR = __spreadArray([], chemComponentsCurrentRep, true);
            var index = chemComponents.indexOf(id);
            compsCR[index] = value;
            setChemComponentsCurrentRep(compsCR);
        }
    };
    var anchorRef = React.useRef();
    var chemRef = React.useRef();
    if (props.aggregate) {
        return (React.createElement(Box, { className: "ParentChem" },
            props.aggregate && (React.createElement(Box, { paddingLeft: 2, paddingRight: 2 },
                React.createElement(Tooltip, { title: "Summary Settings" },
                    React.createElement(Button, { style: { color: "gray" }, ref: anchorRef, onClick: function () { return setSettingsOpen(true); } },
                        React.createElement(SettingsIcon, null),
                        "\u00A0 Settings")),
                React.createElement(SettingsPopover, { open: settingsOpen, setOpen: setSettingsOpen, anchorEl: anchorRef.current, refreshRepList: function () {
                        loadRepList(true);
                    } }),
                React.createElement(Tooltip, { title: "Add Component" },
                    React.createElement(Button, { style: { color: "gray" }, onClick: function () { return addComp(); } },
                        React.createElement(AddCircleOutlineIcon, null),
                        "\u00A0 Add View")))),
            React.createElement("div", { ref: chemRef, className: "chemComponents" },
                chemComponents.length > 1 && (React.createElement("div", { style: {
                        width: (props.rdkitSettings.width + 20) * chemComponents.length,
                    } }, chemComponents.map(function (x, i) {
                    return (React.createElement("div", { key: x, style: {
                            width: props.rdkitSettings.width + 20,
                            float: "left",
                        } },
                        React.createElement(ChemLegend, { chemRef: chemRef, setCurrentRep: function (value) { return setCurrentRep(value, x); }, currentRep: chemComponentsCurrentRep[i], removeComponent: function () { return removeComponent(x); }, id: x, rep_list: repList, selection: props.selection, aggregate: props.aggregate })));
                }))),
                chemComponents.length <= 1 && (React.createElement("div", null,
                    React.createElement("div", { style: { minWidth: props.rdkitSettings.width }, key: chemComponents[0] },
                        React.createElement(ChemLegend, { chemRef: chemRef, setCurrentRep: function (value) {
                                return setCurrentRep(value, chemComponents[0]);
                            }, currentRep: chemComponentsCurrentRep[0], id: chemComponents[0], rep_list: repList, selection: props.selection, aggregate: props.aggregate })))))));
    }
    else {
        return (React.createElement(ChemLegend, { id: -1, rep_list: repList, selection: props.selection, aggregate: props.aggregate }));
    }
});
var loading_area = "chemlegend_loading_area";
var UPDATER = "chemdetail";
var ChemLegend = connector_Chem(/** @class */ (function (_super) {
    __extends(class_1, _super);
    function class_1(props) {
        var _this = _super.call(this, props) || this;
        _this.state = {
            checkedList: [],
        };
        return _this;
    }
    class_1.prototype.render = function () {
        var _this = this;
        var handleMouseEnter = function (i) {
            var hover_item = null;
            if (i >= 0) {
                hover_item = _this.props.selection[i];
            }
            _this.props.setHoverstate(hover_item, UPDATER);
        };
        var handleMouseOut = function () {
            var hover_item = null;
            _this.props.setHoverstate(hover_item, UPDATER);
        };
        var setCheckedList = function (value) {
            var set_val = isFunction(value)
                ? value(_this.state.checkedList)
                : value;
            _this.setState(__assign(__assign({}, _this.state), { checkedList: set_val }));
        };
        var handle_filter = function () {
            var filter_instances = _this.props.selection.filter(function (x, i) { return _this.state.checkedList[i]; });
            if (filter_instances.length > 0) {
                setCheckedList([]);
                _this.props.setCurrentAggregation(filter_instances.map(function (e) { return e.__meta__.meshIndex; }));
            }
            else {
                alert("Please, select at least one Compound in the Summary View to filter.");
            }
        };
        if (this.props.aggregate) {
            return (React.createElement("div", { className: "ParentImg" },
                React.createElement(Box, { paddingLeft: 2, paddingTop: 1, paddingRight: 2 },
                    React.createElement(RepresentationList, { value: this.props.currentRep, onChange: this.props.setCurrentRep, rep_list: this.props.rep_list, hoverSettings: this.props.hoverSettings })),
                React.createElement(Box, { paddingLeft: 2, paddingTop: 1, paddingRight: 2 },
                    React.createElement(Button, { size: "small", variant: "outlined", onClick: function () {
                            handle_filter();
                        } },
                        React.createElement(FilterListIcon, { fontSize: "small" }),
                        "\u00A0Confirm Selection"),
                    this.props.removeComponent && (React.createElement(IconButton, { onClick: this.props.removeComponent },
                        React.createElement(DeleteIcon, null)))),
                React.createElement(LoadingIndicatorView, { area: loading_area + this.props.id }),
                React.createElement(ImageView, { chemRef: this.props.chemRef, id: this.props.id, setCheckedList: setCheckedList, selection: this.props.selection, columns: this.props.columns, aggregate: this.props.aggregate, current_rep: this.props.currentRep, handleMouseEnter: handleMouseEnter, handleMouseOut: handleMouseOut })));
        }
        return (React.createElement("div", null,
            React.createElement(ImageView, { id: this.props.id, selection: this.props.selection, columns: this.props.columns, aggregate: this.props.aggregate })));
    };
    return class_1;
}(React.Component)));
function get_smiles_col(columns) {
    var smiles_col = "SMILES";
    // TODO: find by meta_data -> how to handle multiple smiles columns? for now: just take first column that contains "smiles"
    if (!(smiles_col in columns)) {
        var col_names = Object.keys(columns);
        for (var key in col_names) {
            var col = col_names[key];
            if (col.toLowerCase().includes("smiles")) {
                smiles_col = col;
                break;
            }
        }
    }
    return smiles_col;
}
var mapStateToProps_Img = function (state) { return ({
    datasetState: state.dataset,
    hoverState: state.hoverState,
    rdkitSettings: state.rdkitSettings,
}); };
var mapDispatchToProps_Img = function (dispatch) { return ({}); };
var connector_Img = connect(mapStateToProps_Img, mapDispatchToProps_Img);
function addHighlight(element) {
    if (element && element.style) {
        element.style["border"] = "solid black 1px";
    }
}
function removeHighlight(element) {
    if (element && element.style) {
        element.style["border"] = "solid white 1px";
    }
}
var ImageView = connector_Img(function (_a) {
    var chemRef = _a.chemRef, id = _a.id, hoverState = _a.hoverState, selection = _a.selection, columns = _a.columns, aggregate = _a.aggregate, handleMouseEnter = _a.handleMouseEnter, handleMouseOut = _a.handleMouseOut, current_rep = _a.current_rep, setCheckedList = _a.setCheckedList, rdkitSettings = _a.rdkitSettings, datasetState = _a.datasetState;
    var _b = React.useState(React.createElement("div", null)), comp = _b[0], setComp = _b[1];
    var ref = React.useRef();
    var _c = useCancellablePromise(), cancellablePromise = _c.cancellablePromise, cancelPromises = _c.cancelPromises;
    React.useEffect(function () {
        cancelPromises(); // cancel all unresolved promises
    }, [selection, current_rep]);
    React.useEffect(function () {
        var _a, _b, _c, _d, _e, _f;
        if (setCheckedList)
            setCheckedList([]);
        var smiles_col = get_smiles_col(columns);
        var onUpdateItem = function (i, val) {
            setCheckedList(function (checkedList) {
                var list = checkedList.map(function (item, j) {
                    if (j === i) {
                        return val;
                    }
                    else {
                        return item;
                    }
                });
                return list;
            });
        };
        // TODO: find by meta_data -> how to handle multiple smiles columns?
        // for (const col_name in columns) {
        //     let col = columns[col_name];
        //     if(col.metaInformation.imgSmiles)
        //         smiles_col = col_name;
        // }
        if (smiles_col in columns) {
            setComp(React.createElement("div", null));
            if (selection.length > 0) {
                if (aggregate) {
                    var formData_2 = new FormData();
                    formData_2.append("current_rep", current_rep);
                    selection.forEach(function (row) {
                        formData_2.append("smiles_list", row[smiles_col]);
                    });
                    formData_2.append("contourLines", (_a = rdkitSettings.contourLines) === null || _a === void 0 ? void 0 : _a.toString());
                    formData_2.append("scale", (_b = rdkitSettings.scale) === null || _b === void 0 ? void 0 : _b.toString());
                    formData_2.append("sigma", (_c = rdkitSettings.sigma) === null || _c === void 0 ? void 0 : _c.toString());
                    formData_2.append("showMCS", (_d = rdkitSettings.showMCS) === null || _d === void 0 ? void 0 : _d.toString());
                    formData_2.append("width", (_e = rdkitSettings.width) === null || _e === void 0 ? void 0 : _e.toString());
                    formData_2.append("doAlignment", (_f = rdkitSettings.doAlignment) === null || _f === void 0 ? void 0 : _f.toString());
                    formData_2.append("filename", datasetState.info.path);
                    var controller = new AbortController();
                    trackPromise(cancellablePromise(CIMEBackendFromEnv.getStructuresFromSmilesList(formData_2, controller), controller).then(function (x) {
                        var img_lst = x["img_lst"].map(function (base64, i) {
                            setCheckedList(function (checkedList) {
                                var cpy_checked_list = __spreadArray([], checkedList, true);
                                if (cpy_checked_list.length <= i) {
                                    cpy_checked_list.push(false);
                                }
                                return cpy_checked_list;
                            });
                            return (React.createElement(Grid, { className: "legend_multiple", key: i, item: true },
                                React.createElement(FormControlLabel, { labelPlacement: "bottom", control: React.createElement(Checkbox, { color: "primary", onChange: function (event) {
                                            onUpdateItem(i, event.target.checked);
                                        } }), label: React.createElement("img", { src: "data:image/jpeg;base64," + base64, alt: "", onMouseEnter: function () {
                                            handleMouseEnter(i);
                                        }, onMouseOver: function () {
                                            handleMouseEnter(i);
                                        }, onMouseLeave: function () {
                                            handleMouseOut();
                                        } }) }),
                                React.createElement(Typography, { style: { paddingLeft: 5 }, variant: "subtitle2" },
                                    "ID: ",
                                    selection[i]["ID"])));
                        }); //key={props.selection[i][smiles_col]} --> gives error because sometimes smiles ocure twice
                        //<div dangerouslySetInnerHTML={{ __html: img_lst.join("") }} />
                        setComp(img_lst);
                    }), loading_area + id);
                }
                else {
                    var row_1 = selection[0];
                    var controller = new AbortController();
                    cancellablePromise(CIMEBackendFromEnv.getStructureFromSmiles(datasetState.info.path, row_1[smiles_col], false, controller), controller)
                        .then(function (x) {
                        setComp(React.createElement("div", null,
                            React.createElement("img", { className: "legend_single", src: "data:image/jpeg;base64," + x, alt: "" }),
                            React.createElement(Typography, { style: { paddingLeft: 5 }, variant: "subtitle2" },
                                "ID: ",
                                row_1["ID"])));
                    })
                        .catch(function (error) { return console.log(error); });
                }
            }
            else {
                setComp(React.createElement("div", null, "No Selection"));
            }
        }
        else {
            setComp(React.createElement("div", null, "No SMILES column found"));
        }
    }, [selection]);
    React.useEffect(function () {
        var _a, _b, _c, _d, _e, _f;
        if (aggregate) {
            var smiles_col_2 = get_smiles_col(columns);
            if (smiles_col_2 in columns) {
                var imgList_1 = ref.current.childNodes;
                if (selection.length === imgList_1.length) {
                    ref.current.style.display = "none";
                    var formData_3 = new FormData();
                    formData_3.append("current_rep", current_rep);
                    selection.forEach(function (row) {
                        formData_3.append("smiles_list", row[smiles_col_2]);
                    });
                    formData_3.append("contourLines", (_a = rdkitSettings.contourLines) === null || _a === void 0 ? void 0 : _a.toString());
                    formData_3.append("scale", (_b = rdkitSettings.scale) === null || _b === void 0 ? void 0 : _b.toString());
                    formData_3.append("sigma", (_c = rdkitSettings.sigma) === null || _c === void 0 ? void 0 : _c.toString());
                    formData_3.append("showMCS", (_d = rdkitSettings.showMCS) === null || _d === void 0 ? void 0 : _d.toString());
                    formData_3.append("width", (_e = rdkitSettings.width) === null || _e === void 0 ? void 0 : _e.toString());
                    formData_3.append("doAlignment", (_f = rdkitSettings.doAlignment) === null || _f === void 0 ? void 0 : _f.toString());
                    formData_3.append("filename", datasetState.info.path);
                    var controller = new AbortController();
                    trackPromise(cancellablePromise(CIMEBackendFromEnv.getStructuresFromSmilesList(formData_3, controller), controller).then(function (x) {
                        x["img_lst"].forEach(function (base64, i) {
                            var cur_img = imgList_1[i].getElementsByTagName("img")[0];
                            cur_img.src = "data:image/jpeg;base64," + base64;
                            // cur_img.width = props.rdkitSettings.width;
                            // cur_img.height = props.rdkitSettings.width;
                        });
                        ref.current.style.display = "flex";
                    }), loading_area + id);
                }
            }
        }
    }, [current_rep, rdkitSettings.refresh]);
    React.useEffect(function () {
        if (aggregate) {
            //@ts-ignore
            var container = chemRef === null || chemRef === void 0 ? void 0 : chemRef.current;
            var imgContainer = container.getElementsByClassName("chem-grid")[0];
            //@ts-ignore
            var imgList = imgContainer.childNodes;
            if (hoverState && hoverState.data) {
                var idx = selection.findIndex(function (x) {
                    return x &&
                        x["__meta__"] &&
                        hoverState.data["__meta__"] &&
                        x["__meta__"]["meshIndex"] ==
                            hoverState.data["__meta__"]["meshIndex"];
                });
                if (idx >= 0 && imgList.length > 0) {
                    for (var i in imgList) {
                        var img_div = imgList[i];
                        removeHighlight(img_div);
                    }
                    addHighlight(imgList[idx]);
                    if (hoverState.updater != UPDATER) {
                        if (container && imgList[idx]) {
                            //@ts-ignore
                            container.scrollTop =
                                imgList[idx].offsetTop - container.offsetTop;
                            // imgList[idx]?.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'start' }) // this seems to be buggy sometimes
                        }
                    }
                }
            }
            else {
                for (var i in imgList) {
                    var img_div = imgList[i];
                    removeHighlight(img_div);
                }
            }
        }
    }, [hoverState.data, hoverState.updater]);
    return (React.createElement("div", { className: "chemContainer" },
        React.createElement(Grid, { ref: ref, className: "chem-grid", container: true }, comp)));
});
function ValueLabelComponent(props) {
    var children = props.children, open = props.open, value = props.value;
    return (React.createElement(Tooltip, { open: open, enterTouchDelay: 0, placement: "top", title: value }, children));
}
var mapStateToProps_settings = function (state) { return ({
    rdkitSettings: state.rdkitSettings,
}); };
var mapDispatchToProps_settings = function (dispatch) { return ({
    setContourLines: function (input) { return dispatch(setRDKit_contourLines(input)); },
    setScale: function (input) { return dispatch(setRDKit_scale(input)); },
    setSigma: function (input) { return dispatch(setRDKit_sigma(input)); },
    setShowMCS: function (input) { return dispatch(setRDKit_showMCS(input)); },
    setWidth: function (input) { return dispatch(setRDKit_width(input)); },
    setRefresh: function (input) { return dispatch(setRDKit_refresh(input)); },
    setDoAlignment: function (input) { return dispatch(setRDKit_doAlignment(input)); },
}); };
var connector_settings = connect(mapStateToProps_settings, mapDispatchToProps_settings);
var SettingsPopover = connector_settings(function (_a) {
    var open = _a.open, setOpen = _a.setOpen, anchorEl = _a.anchorEl, refreshRepList = _a.refreshRepList, rdkitSettings = _a.rdkitSettings, setContourLines = _a.setContourLines, setScale = _a.setScale, setSigma = _a.setSigma, setShowMCS = _a.setShowMCS, setWidth = _a.setWidth, setRefresh = _a.setRefresh, setDoAlignment = _a.setDoAlignment;
    return (React.createElement(Popover, { disablePortal: true, id: "dialog to open", open: open, anchorEl: anchorEl, onClose: function () { return setOpen(function () { return false; }); }, anchorOrigin: {
            vertical: "bottom",
            horizontal: "right",
        }, transformOrigin: {
            vertical: "bottom",
            horizontal: "left",
        } },
        React.createElement("div", null,
            React.createElement(Paper, { style: { padding: 10, minWidth: 300 } },
                React.createElement(Grid, { container: true, spacing: 2, direction: 'column', justifyContent: 'stretch', sx: {
                        '& .MuiTextField-root': { width: '200px' },
                    } },
                    React.createElement(Grid, { item: true },
                        React.createElement(Button, { size: "small", variant: "outlined", "aria-label": "Refresh Representation List", onClick: function () { return refreshRepList(true); } },
                            React.createElement(RefreshIcon, null),
                            "Refresh Representation List")),
                    React.createElement(Grid, { item: true },
                        React.createElement(Typography, { variant: "subtitle2", gutterBottom: true }, "RDKit Settings")),
                    React.createElement(Grid, { item: true, container: true, alignItems: 'center' },
                        React.createElement(TextField, { label: "Contour Lines", type: "number", InputLabelProps: {
                                shrink: true
                            }, value: rdkitSettings.contourLines, onChange: function (event) {
                                var val = parseInt(event.target.value);
                                if (isNaN(val))
                                    setContourLines(event.target.value);
                                else
                                    setContourLines(Math.max(val, 0));
                            } }),
                        React.createElement("div", { style: { flexGrow: 1 } }),
                        React.createElement(Tooltip, { title: "Number of Contour Lines [0; \u221E] \u2208 \u2115" },
                            React.createElement(InfoIcon, { fontSize: "small" }))),
                    React.createElement(Grid, { item: true, container: true, alignItems: 'center' },
                        React.createElement(TextField, { label: "Scale", type: "number", InputLabelProps: {
                                shrink: true
                            }, value: rdkitSettings.scale, onChange: function (event) {
                                var val = parseFloat(event.target.value);
                                if (isNaN(val))
                                    setScale(event.target.value);
                                else
                                    setScale(Math.max(val, -1));
                            } }),
                        React.createElement("div", { style: { flexGrow: 1 } }),
                        React.createElement(Tooltip, { title: "Weight Scale [-1; \u221E] \u2208 \u211D" },
                            React.createElement(InfoIcon, { fontSize: "small" }))),
                    React.createElement(Grid, { item: true, container: true, alignItems: 'center' },
                        React.createElement(TextField, { label: "Sigma", type: "number", InputLabelProps: {
                                shrink: true
                            }, inputProps: {
                                step: 0.1
                            }, value: rdkitSettings.sigma, onChange: function (event) {
                                var val = parseFloat(event.target.value);
                                if (isNaN(val))
                                    setSigma(event.target.value);
                                else
                                    setSigma(Math.max(val, 0));
                            } }),
                        React.createElement("div", { style: { flexGrow: 1 } }),
                        React.createElement(Tooltip, { title: "Sigma for Gaussian ]0; \u221E] \u2208 \u211D. Default of 0 signals the algorithm to infer the value." },
                            React.createElement(InfoIcon, { fontSize: "small" }))),
                    React.createElement(Grid, { item: true },
                        React.createElement(FormControlLabel, { control: React.createElement(Switch, { color: "primary", checked: rdkitSettings.showMCS, onChange: function (_, value) {
                                    setShowMCS(value);
                                } }), label: "Show MCS" })),
                    React.createElement(Grid, { item: true },
                        React.createElement(FormControlLabel, { control: React.createElement(Switch, { color: "primary", checked: rdkitSettings.doAlignment, onChange: function (_, value) {
                                    setDoAlignment(value);
                                } }), label: "Align Structure" })),
                    React.createElement(Grid, { item: true },
                        React.createElement(Typography, { style: { paddingTop: 10 }, gutterBottom: true }, "Image Width"),
                        React.createElement(Slider, { value: rdkitSettings.width, onChange: function (event, new_val) {
                                setWidth(new_val);
                            }, min: 50, max: 500, step: 10 })),
                    React.createElement(Grid, { item: true },
                        React.createElement(Button, { style: { marginTop: 3, maxWidth: 150 }, size: "small", variant: "outlined", onClick: function () {
                                setRefresh((rdkitSettings.refresh += 1));
                            } }, "Apply Settings")))))));
});
var RepresentationList = function (props) {
    var options = props.rep_list.map(function (rep) {
        var split = rep.split("_");
        var inputVal = split.pop();
        var group = split.join("_");
        group = group.replace("atom.dprop.", "");
        group = group.replace("atom.dprop", "");
        return {
            group: group,
            value: rep,
            inputValue: inputVal,
        };
    });
    var filterOptions = createFilterOptions({
        stringify: function (option) {
            return option.value;
        },
    });
    return (React.createElement(Autocomplete, { size: "small", className: props.className, filterOptions: filterOptions, onChange: function (event, newValue) {
            if (newValue)
                props.onChange(newValue.value);
        }, disablePortal: props.hoverSettings.windowMode == WindowMode.Extern, options: options.sort(function (a, b) { return -b.group.localeCompare(a.group); }), groupBy: function (option) { return option.group; }, getOptionLabel: function (option) { return option.inputValue; }, style: { maxWidth: 300 }, defaultValue: options[0], renderInput: function (params) { return (React.createElement(TextField, __assign({}, params, { label: "Choose Representation", variant: "outlined" }))); }, isOptionEqualToValue: function (option, value) { return option.value == value.value; } }));
};
//getOptionSelected={(option:any, value) => {return option.value == value.value;}}
//# sourceMappingURL=ChemDetail.js.map