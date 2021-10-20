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
import { useState } from "react";
import { PSEContextProvider, API, Application, PluginRegistry, createRootReducer, } from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { CIMEReducers } from "./State/Store";
export var DEMO = false;
PluginRegistry.getInstance().registerPlugin(new ChemPlugin());
export function CIMEApp() {
    var context = useState(new API(null, createRootReducer(CIMEReducers)))[0];
    return (_jsx(PSEContextProvider, __assign({ context: context }, { children: _jsx(Application, { config: {}, features: {
                disableEmbeddings: {
                    tsne: true,
                    forceatlas: true,
                },
            }, overrideComponents: {
                datasetTab: DatasetTabPanel,
                appBar: CimeAppBar,
                tabs: [
                    {
                        name: "lineup",
                        //@ts-ignore
                        tab: LineUpTabPanel,
                        title: "LineUp Integration",
                        description: "Settings for LineUp Integration",
                        icon: null,
                    },
                ],
                detailViews: [
                    {
                        name: "lineup",
                        //@ts-ignore
                        view: LineUpContext,
                    },
                ],
            } }, void 0) }), void 0));
}
//# sourceMappingURL=CIMEApp.js.map