var __rest = (this && this.__rest) || function (s, e) {
    var t = {};
    for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p) && e.indexOf(p) < 0)
        t[p] = s[p];
    if (s != null && typeof Object.getOwnPropertySymbols === "function")
        for (var i = 0, p = Object.getOwnPropertySymbols(s); i < p.length; i++) {
            if (e.indexOf(p[i]) < 0 && Object.prototype.propertyIsEnumerable.call(s, p[i]))
                t[p[i]] = s[p[i]];
        }
    return t;
};
import * as React from "react";
import { PSEContextProvider, API, Application, PluginRegistry, createRootReducer, DEFAULT_UMAP_SETTINGS, } from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { CIMEReducers } from "./State/Store";
import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
import "./LineUpContext.scss";
export var DEMO = false;
PluginRegistry.getInstance().registerPlugin(new ChemPlugin());
export var DEFAULT_CIME_APP_CONFIG = {
    config: {
        preselect: {
            initOnMount: false
        }
    },
    features: {
        embeddings: [
            { id: "umap", name: "UMAP", settings: DEFAULT_UMAP_SETTINGS },
        ],
        showSummaryAttributes: false
    },
    overrideComponents: {
        datasetTab: DatasetTabPanel,
        appBar: CimeAppBar,
    }
};
// TODO: @moritz We are currently using the global object directly, ideally we make it passable as prop.
export var CIMEAppContext = new API(undefined, createRootReducer(CIMEReducers));
export function CIMEApp(_a) {
    var _b = _a.providePSEContext, providePSEContext = _b === void 0 ? true : _b, props = __rest(_a, ["providePSEContext"]);
    var _c = React.useState(null), merged = _c[0], setMerged = _c[1];
    React.useEffect(function () {
        setMerged(merge(cloneDeep(DEFAULT_CIME_APP_CONFIG), props));
    }, [props.config, props.features, props.overrideComponents]);
    console.log(merged);
    if (!merged) {
        return null;
    }
    var component = React.createElement("div", { style: { width: '100%', height: '100%' } },
        React.createElement(Application, { ref: props.pseRef, config: merged === null || merged === void 0 ? void 0 : merged.config, features: merged === null || merged === void 0 ? void 0 : merged.features, 
            //@ts-ignore
            overrideComponents: merged === null || merged === void 0 ? void 0 : merged.overrideComponents }));
    return (providePSEContext ?
        React.createElement(PSEContextProvider, { context: CIMEAppContext }, component) : component);
}
//# sourceMappingURL=CIMEApp.js.map