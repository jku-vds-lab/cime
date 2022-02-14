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
    var renderMenu = (React.createElement(Menu, { anchorEl: anchorEl, anchorOrigin: {
            vertical: 'top',
            horizontal: 'right',
        }, id: menuId, keepMounted: true, transformOrigin: {
            vertical: 'top',
            horizontal: 'right',
        }, open: isMenuOpen, onClose: handleMenuClose },
        React.createElement(MenuItem, { onClick: handleMenuClose }, "About"),
        React.createElement(MenuItem, { onClick: handleMenuClose }, "Help"),
        React.createElement(MenuItem, { onClick: handleMenuClose }, "...")));
    return (React.createElement(Box, { sx: { flexGrow: 1 } },
        React.createElement(AppBar, { variant: "outlined", elevation: 0, position: "relative", color: "transparent" },
            React.createElement(Toolbar, null,
                React.createElement("a", { href: "https://www.bayer.com", target: "_blank", rel: "noreferrer" },
                    React.createElement("img", { style: { height: 48 }, src: BayerLogo, alt: "Powered By Bayer" })),
                React.createElement(Typography, { variant: "h6", style: { marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" } }, "CIME: ChemInformatics Model Explorer"),
                React.createElement(Box, { sx: { flexGrow: 1 } }),
                React.createElement(IconButton, { size: "large", edge: "end", "aria-label": "account of current user", "aria-controls": menuId, onClick: handleProfileMenuOpen, "aria-haspopup": "true", color: "inherit" },
                    React.createElement(MoreIcon, null)))),
        renderMenu));
}
//# sourceMappingURL=CimeAppBar.js.map