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
import * as React from "react";
import { Box, IconButton, Menu, MenuItem, Toolbar } from "@mui/material";
import { AppBar, Typography } from "@mui/material";
import BayerLogo from "../assets/bayer-logo.svg";
import MoreIcon from '@mui/icons-material/MoreVert';
export function CimeAppBar() {
    var menuId = 'primary-search-account-menu';
    var _a = React.useState(null), anchorEl = _a[0], setAnchorEl = _a[1];
    var isMenuOpen = Boolean(anchorEl);
    var handleMenuClose = function () {
        setAnchorEl(null);
    };
    var handleProfileMenuOpen = function (event) {
        setAnchorEl(event.currentTarget);
    };
    var renderMenu = (_jsxs(Menu, __assign({ anchorEl: anchorEl, anchorOrigin: {
            vertical: 'top',
            horizontal: 'right',
        }, id: menuId, keepMounted: true, transformOrigin: {
            vertical: 'top',
            horizontal: 'right',
        }, open: isMenuOpen, onClose: handleMenuClose }, { children: [_jsx(MenuItem, __assign({ onClick: handleMenuClose }, { children: "About" }), void 0), _jsx(MenuItem, __assign({ onClick: handleMenuClose }, { children: "Help" }), void 0), _jsx(MenuItem, __assign({ onClick: handleMenuClose }, { children: "..." }), void 0)] }), void 0));
    return (_jsxs(Box, __assign({ sx: { flexGrow: 1 } }, { children: [_jsx(AppBar, __assign({ variant: "outlined", elevation: 0, position: "relative", color: "transparent" }, { children: _jsxs(Toolbar, { children: [_jsx("a", __assign({ href: "https://www.bayer.com", target: "_blank", rel: "noreferrer" }, { children: _jsx("img", { style: { height: 48 }, src: BayerLogo, alt: "Powered By Bayer" }, void 0) }), void 0), _jsx(Typography, __assign({ variant: "h6", style: { marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" } }, { children: "CIME: ChemInformatics Model Explorer" }), void 0), _jsx(Box, { sx: { flexGrow: 1 } }, void 0), _jsx(IconButton, __assign({ size: "large", edge: "end", "aria-label": "account of current user", "aria-controls": menuId, onClick: handleProfileMenuOpen, "aria-haspopup": "true", color: "inherit" }, { children: _jsx(MoreIcon, {}, void 0) }), void 0)] }, void 0) }), void 0), renderMenu] }), void 0));
}
//# sourceMappingURL=CimeAppBar.js.map