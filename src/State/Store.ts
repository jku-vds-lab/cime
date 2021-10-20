import { RootState } from "projection-space-explorer";
import { combineReducers } from "redux";
import lineUpInput from "./LineUpInputDuck";
import rdkitSettings from "./RDKitSettingsDuck";

export const CIMEReducers = {
  lineUpInput: lineUpInput,
  rdkitSettings: rdkitSettings,
};

const combined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof combined>;

export type AppState = RootState & CimeState;
