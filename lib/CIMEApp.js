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
import { jsx as _jsx } from "react/jsx-runtime";
import * as React from "react";
import { PSEContextProvider, API, Application, PluginRegistry, createRootReducer, PSEIcons, DEFAULT_UMAP_SETTINGS, } from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { CIMEReducers } from "./State/Store";
import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
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
    },
    overrideComponents: {
        datasetTab: DatasetTabPanel,
        appBar: CimeAppBar,
        tabs: [
            {
                name: "lineup",
                // @ts-ignore TODO: @moritz error after I added the correct typing for DEFAULT_CIME_APP_CONFIG
                tab: LineUpTabPanel,
                title: "LineUp Integration",
                description: "Settings for LineUp Integration",
                icon: PSEIcons.PseLineup,
            },
        ],
        detailViews: [
            {
                name: "lineup",
                // @ts-ignore TODO: @moritz error after I added the correct typing for DEFAULT_CIME_APP_CONFIG
                view: LineUpContext,
            },
        ],
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
    var component = _jsx("div", __assign({ style: { width: '100%', height: '100%' } }, { children: _jsx(Application, { ref: props.pseRef, config: merged === null || merged === void 0 ? void 0 : merged.config, features: merged === null || merged === void 0 ? void 0 : merged.features, 
            //@ts-ignore
            overrideComponents: merged === null || merged === void 0 ? void 0 : merged.overrideComponents }, void 0) }), void 0);
    return (providePSEContext ?
        _jsx(PSEContextProvider, __assign({ context: CIMEAppContext }, { children: component }), void 0) : component);
}
//# sourceMappingURL=CIMEApp.js.map