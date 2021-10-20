/**
 * Duck file for the LineUp input data
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
import { filter } from "lodash";
// const SET_DATA = "ducks/lineUpInput/SET_DATA"
// const SET_COLUMNS = "ducks/lineUpInput/SET_COLUMNS"
var SET_VISIBILITY = "ducks/lineUpInput/SET_VISIBILITY";
var SET_DUMP = "ducks/lineUpInput/SET_DUMP";
var SET_FILTER = "ducks/lineUpInput/SET_FILTER";
var UPDATE_FILTER = "ducks/lineUpInput/UPDATE_FILTER";
var SET_LINEUP = "ducks/lineUpInput/SET_LINEUP";
var SET_UPDATE = "ducks/lineUpInput/SET_UPDATE";
// export const setLineUpInput_data = input => ({
//     type: SET_DATA,
//     input: input
// });
// export const setLineUpInput_columns = input => ({
//     type: SET_COLUMNS,
//     input: input
// });
export var setLineUpInput_visibility = function (input) { return ({
    type: SET_VISIBILITY,
    input: input,
}); };
export var setLineUpInput_dump = function (input) { return ({
    type: SET_DUMP,
    input: input,
}); };
export var setLineUpInput_filter = function (input) { return ({
    type: SET_FILTER,
    input: input,
}); };
export var updateLineUpInput_filter = function (input) { return ({
    type: UPDATE_FILTER,
    input: input,
}); };
export var setLineUpInput_lineup = function (input) { return ({
    type: SET_LINEUP,
    input: input,
}); };
export var setLineUpInput_update = function (input) { return ({
    type: SET_UPDATE,
    input: input,
}); };
var initialState = {
    // data: null,
    // columns: null,
    show: false,
    dump: "",
    filter: null,
    previousfilter: null,
    lineup: null,
    update: 0,
};
var lineUpInput = function (state, action) {
    if (state === void 0) { state = initialState; }
    switch (action.type) {
        // case SET_DATA:
        //     return {...state, data: action.input}
        // case SET_COLUMNS:
        //     return {...state, columns: action.input}
        case SET_VISIBILITY:
            return __assign(__assign({}, state), { show: action.input });
        case SET_DUMP:
            return __assign(__assign({}, state), { dump: action.input });
        case SET_FILTER:
            var prev_filter = __assign({}, state.filter);
            return __assign(__assign({}, state), { previousfilter: prev_filter, filter: action.input });
        case UPDATE_FILTER:
            if (state.filter &&
                Object.keys(state.filter).includes(action.input["key"])) {
                if (state.filter[action.input["key"]] == action.input["val_old"]) {
                    var filter_new = __assign({}, filter);
                    filter_new[action.input["key"]] = action.input["val_new"];
                    return __assign(__assign({}, state), { filter: filter_new });
                }
            }
            return state;
        case SET_LINEUP:
            return __assign(__assign({}, state), { lineup: action.input });
        case SET_UPDATE:
            var cur = state.update;
            return __assign(__assign({}, state), { update: cur + 1 });
        default:
            return state;
    }
};
export default lineUpInput;
//# sourceMappingURL=LineUpInputDuck.js.map