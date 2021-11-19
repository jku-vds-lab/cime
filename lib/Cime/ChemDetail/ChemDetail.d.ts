/// <reference types="react" />
import "./chem.scss";
import { ConnectedProps } from "react-redux";
declare const connector_Chem: import("react-redux").InferableComponentEnhancerWithProps<{
    dataset: import("projection-space-explorer").Dataset;
    hoverSettings: {
        windowMode: any;
    };
    rdkitSettings: import("../../State/RDKitSettingsDuck").RDKitSettingsType;
    columns: {
        [name: string]: {
            distinct: any;
            isNumeric: boolean;
            metaInformation: any;
            featureType: import("projection-space-explorer").FeatureType;
            range: any;
            featureLabel: string;
            project: boolean;
        };
    };
} & {
    setCurrentAggregation: (samples: number[]) => any;
    setHoverstate: (state: any, updater: any) => any;
}, {}>;
/**
 * Type that holds the props we declared above in mapStateToProps and mapDispatchToProps
 */
declare type PropsFromRedux_Chem = ConnectedProps<typeof connector_Chem>;
declare type Props_Chem_Parent = PropsFromRedux_Chem & {
    selection: any[];
    aggregate: boolean;
    mcs_only?: boolean;
    diff?: boolean;
    selection_ref?: any[];
};
export declare const ChemLegendParent: import("react-redux").ConnectedComponent<(props: Props_Chem_Parent) => JSX.Element, import("react-redux").Omit<{
    dataset: import("projection-space-explorer").Dataset;
    hoverSettings: {
        windowMode: any;
    };
    rdkitSettings: import("../../State/RDKitSettingsDuck").RDKitSettingsType;
    columns: {
        [name: string]: {
            distinct: any;
            isNumeric: boolean;
            metaInformation: any;
            featureType: import("projection-space-explorer").FeatureType;
            range: any;
            featureLabel: string;
            project: boolean;
        };
    };
} & {
    setCurrentAggregation: (samples: number[]) => any;
    setHoverstate: (state: any, updater: any) => any;
} & {
    selection: any[];
    aggregate: boolean;
    mcs_only?: boolean;
    diff?: boolean;
    selection_ref?: any[];
}, "columns" | "rdkitSettings" | "dataset" | "hoverSettings" | "setCurrentAggregation" | "setHoverstate">>;
export {};
//# sourceMappingURL=ChemDetail.d.ts.map