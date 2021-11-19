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
import { jsx as _jsx } from "react/jsx-runtime";
import * as React from "react";
import { PSEContextProvider, API, Application, PluginRegistry, createRootReducer, PSEIcons, } from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { CIMEReducers } from "./State/Store";
import merge from 'lodash/merge';
import clone from "lodash/clone";
export var DEMO = false;
PluginRegistry.getInstance().registerPlugin(new ChemPlugin());
export var DEFAULT_CIME_APP_CONFIG = {
    config: {
        preselect: {
            initOnMount: false
        }
    },
    features: {
        disableEmbeddings: {
            tsne: true,
            forceatlas: true,
        },
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
export function CIMEApp(props) {
    var _a = React.useState(null), merged = _a[0], setMerged = _a[1];
    React.useEffect(function () {
        setMerged(merge(clone(DEFAULT_CIME_APP_CONFIG), props));
    }, [props.config, props.features, props.overrideComponents]);
    return (merged ? _jsx(PSEContextProvider, __assign({ context: CIMEAppContext }, { children: _jsx("div", __assign({ style: { width: '100%', height: '100%' } }, { children: _jsx(Application, { ref: props.pseRef, config: merged.config, features: merged.features, 
                //@ts-ignore
                overrideComponents: merged.overrideComponents }, void 0) }), void 0) }), void 0) : null);
}
//# sourceMappingURL=CIMEApp.js.map