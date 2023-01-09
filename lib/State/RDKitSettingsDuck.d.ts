/**
 * Duck file for the RDKit Settings
 */
export declare const setRDKit_contourLines: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_scale: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_sigma: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_refresh: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_showMCS: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_width: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_doAlignment: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_colorDomain: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_showBarChart: (input: any) => {
    type: string;
    input: any;
};
export declare const setRDKit_showAttributes: (input: any) => {
    type: string;
    input: any;
};
export type RDKitSettingsType = {
    contourLines: number;
    scale: number;
    sigma: number;
    refresh: number;
    showMCS: boolean;
    width: number;
    doAlignment: boolean;
    domain: {
        type: number;
        value: number[];
        deadzone: number;
        thresholds: number[];
    };
    showBarChart: boolean;
    showAttributes: boolean;
};
declare const rdkitSettings: (state: RDKitSettingsType, action: any) => RDKitSettingsType;
export default rdkitSettings;
//# sourceMappingURL=RDKitSettingsDuck.d.ts.map