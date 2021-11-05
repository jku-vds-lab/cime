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
import { useState } from "react";
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
                tab: LineUpTabPanel,
                title: "LineUp Integration",
                description: "Settings for LineUp Integration",
                icon: PSEIcons.PseLineup,
            },
        ],
        detailViews: [
            {
                name: "lineup",
                view: LineUpContext,
            },
        ],
    }
};
export function CIMEApp(props) {
    var context = useState(new API(null, createRootReducer(CIMEReducers)))[0];
    var _a = React.useState(merge(clone(DEFAULT_CIME_APP_CONFIG), props)), merged = _a[0], setMerged = _a[1];
    React.useEffect(function () {
        var cimeAppConfig = clone(DEFAULT_CIME_APP_CONFIG);
        setMerged(merge(cimeAppConfig, props));
    }, [props.config, props.features, props.overrideComponents]);
    return (_jsx(PSEContextProvider, __assign({ context: context }, { children: _jsx(Application, { config: merged.config, features: merged.features, 
            //@ts-ignore
            overrideComponents: merged.overrideComponents }, void 0) }), void 0));
}
//# sourceMappingURL=CIMEApp.js.map