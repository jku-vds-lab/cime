/**
 * Duck file for the RDKit Settings
 */
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
var SET_CONTOURLINES = "ducks/rdkitsettings/SET_CONTOURLINES";
var SET_SCALE = "ducks/rdkitsettings/SET_SCALE";
var SET_SIGMA = "ducks/rdkitsettings/SET_SIGMA";
var SET_REFRESH = "ducks/rdkitsettings/SET_REFRESH";
var SET_SHOW_MCS = "ducks/rdkitsettings/SET_SHOW_MCS";
var SET_WIDTH = "ducks/rdkitsettings/SET_WIDTH";
var SET_DO_ALIGNMENT = "ducks/rdkitsettings/SET_DO_ALIGNMENT";
var SET_DO_DOMAIN = "ducks/rdkitsettings/SET_DO_DOMAIN";
export var setRDKit_contourLines = function (input) { return ({
    type: SET_CONTOURLINES,
    input: input,
}); };
export var setRDKit_scale = function (input) { return ({
    type: SET_SCALE,
    input: input,
}); };
export var setRDKit_sigma = function (input) { return ({
    type: SET_SIGMA,
    input: input,
}); };
export var setRDKit_refresh = function (input) { return ({
    type: SET_REFRESH,
    input: input,
}); };
export var setRDKit_showMCS = function (input) { return ({
    type: SET_SHOW_MCS,
    input: input,
}); };
export var setRDKit_width = function (input) { return ({
    type: SET_WIDTH,
    input: input,
}); };
export var setRDKit_doAlignment = function (input) { return ({
    type: SET_DO_ALIGNMENT,
    input: input,
}); };
export var setRDKit_colorDomain = function (input) { return ({
    type: SET_DO_DOMAIN,
    input: input,
}); };
var initialState = {
    contourLines: 10,
    scale: -1,
    sigma: 0,
    refresh: 0,
    showMCS: false,
    width: 250,
    doAlignment: true,
    domain: {
        type: 0,
        value: null,
        deadzone: 0,
        thresholds: [],
    },
};
var rdkitSettings = function (state, action) {
    if (state === void 0) { state = initialState; }
    switch (action.type) {
        case SET_CONTOURLINES:
            return __assign(__assign({}, state), { contourLines: action.input });
        case SET_SCALE:
            return __assign(__assign({}, state), { scale: action.input });
        case SET_SIGMA:
            return __assign(__assign({}, state), { sigma: action.input });
        case SET_REFRESH:
            return __assign(__assign({}, state), { refresh: action.input });
        case SET_SHOW_MCS:
            return __assign(__assign({}, state), { showMCS: action.input });
        case SET_WIDTH:
            return __assign(__assign({}, state), { width: action.input });
        case SET_DO_ALIGNMENT:
            return __assign(__assign({}, state), { doAlignment: action.input });
        case SET_DO_DOMAIN:
            return __assign(__assign({}, state), { domain: __assign(__assign({}, state.domain), action.input) });
        default:
            return state;
    }
};
export default rdkitSettings;
//# sourceMappingURL=RDKitSettingsDuck.js.map