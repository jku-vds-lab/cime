import { RootState } from "projection-space-explorer";
export declare const CIMEReducers: {
    lineUpInput: (state: import("./LineUpInputDuck").LineUpType, action: any) => import("./LineUpInputDuck").LineUpType;
    rdkitSettings: (state: import("./RDKitSettingsDuck").RDKitSettingsType, action: any) => import("./RDKitSettingsDuck").RDKitSettingsType;
};
declare const combined: import("redux").Reducer<import("redux").CombinedState<{
    lineUpInput: import("./LineUpInputDuck").LineUpType;
    rdkitSettings: import("./RDKitSettingsDuck").RDKitSettingsType;
}>, import("redux").AnyAction>;
/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof combined>;
export type AppState = RootState & CimeState;
export {};
//# sourceMappingURL=Store.d.ts.map