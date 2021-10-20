import { combineReducers } from "redux";
import lineUpInput from "./LineUpInputDuck";
import rdkitSettings from "./RDKitSettingsDuck";
export var CIMEReducers = {
    lineUpInput: lineUpInput,
    rdkitSettings: rdkitSettings,
};
var combined = combineReducers(CIMEReducers);
//# sourceMappingURL=Store.js.map