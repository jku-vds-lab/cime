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
import { Toolbar } from "@mui/material";
import { AppBar, Typography } from "@mui/material";
import * as BayerLogo from "../assets/bayer-logo.svg";
export function CimeAppBar() {
    return (_jsx(AppBar, __assign({ variant: "outlined", elevation: 0, position: "relative", color: "transparent" }, { children: _jsxs(Toolbar, { children: [_jsx("a", __assign({ href: "https://www.bayer.com", target: "_blank", rel: "noreferrer" }, { children: _jsx("img", { style: { height: 48 }, src: BayerLogo, alt: "Powered By Bayer" }, void 0) }), void 0), _jsx(Typography, __assign({ variant: "h6", style: { marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" } }, { children: "CIME: Chem-Informatics Model Explorer" }), void 0)] }, void 0) }), void 0));
}
//# sourceMappingURL=CimeAppBar.js.map