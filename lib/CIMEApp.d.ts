import { BaseConfig, FeatureConfig, ComponentConfig } from "projection-space-explorer";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
export declare const DEMO = false;
export declare type CIMEAppProps = {
    config?: BaseConfig;
    features?: FeatureConfig;
    overrideComponents?: ComponentConfig;
};
export declare const DEFAULT_CIME_APP_CONFIG: {
    config: {
        preselect: {
            initOnMount: boolean;
        };
    };
    features: {
        disableEmbeddings: {
            tsne: boolean;
            forceatlas: boolean;
        };
    };
    overrideComponents: {
        datasetTab: typeof DatasetTabPanel;
        appBar: typeof CimeAppBar;
        tabs: {
            name: string;
            tab: import("react-redux").ConnectedComponent<({ setLineUpInput_visibility, setLineUpInput_filter, lineUpInput, dataset, currentAggregation, splitRef, }: {
                dataset: import("projection-space-explorer").Dataset;
                currentAggregation: {
                    aggregation: number[];
                    selectedClusters: string[];
                    source: "sample" | "cluster";
                };
                lineUpInput: import("./State/LineUpInputDuck").LineUpType;
            } & {
                setLineUpInput_visibility: (value: any) => any;
                setLineUpInput_filter: (value: any) => any;
            } & {
                splitRef: any;
            }) => JSX.Element, import("react-redux").Omit<{
                dataset: import("projection-space-explorer").Dataset;
                currentAggregation: {
                    aggregation: number[];
                    selectedClusters: string[];
                    source: "sample" | "cluster";
                };
                lineUpInput: import("./State/LineUpInputDuck").LineUpType;
            } & {
                setLineUpInput_visibility: (value: any) => any;
                setLineUpInput_filter: (value: any) => any;
            } & {
                splitRef: any;
            }, "lineUpInput" | "dataset" | "currentAggregation" | "setLineUpInput_visibility" | "setLineUpInput_filter">>;
            title: string;
            description: string;
            icon: any;
        }[];
        detailViews: {
            name: string;
            view: import("react-redux").ConnectedComponent<({ dataset, lineUpInput, lineUpInput_data, lineUpInput_columns, currentAggregation, channelColor, setCurrentAggregation, setLineUpInput_lineup, setLineUpInput_visibility, onFilter, activeStory, pointColorScale, setHoverstate, detailView, }: {
                dataset: import("projection-space-explorer").Dataset;
                lineUpInput: import("./State/LineUpInputDuck").LineUpType;
                lineUpInput_data: import("projection-space-explorer").IVector[];
                lineUpInput_columns: {
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
                currentAggregation: {
                    aggregation: number[];
                    selectedClusters: string[];
                    source: "sample" | "cluster";
                };
                activeStory: import("projection-space-explorer").IBook;
                pointColorScale: any;
                channelColor: any;
                detailView: {
                    open: boolean;
                    active: string;
                };
            } & {
                setCurrentAggregation: (samples: number[]) => any;
                setLineUpInput_visibility: (visibility: any) => any;
                setLineUpInput_lineup: (input: any) => any;
                setHoverstate: (state: any, updater: any) => any;
            } & {
                onFilter: any;
            }) => JSX.Element, import("react-redux").Omit<{
                dataset: import("projection-space-explorer").Dataset;
                lineUpInput: import("./State/LineUpInputDuck").LineUpType;
                lineUpInput_data: import("projection-space-explorer").IVector[];
                lineUpInput_columns: {
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
                currentAggregation: {
                    aggregation: number[];
                    selectedClusters: string[];
                    source: "sample" | "cluster";
                };
                activeStory: import("projection-space-explorer").IBook;
                pointColorScale: any;
                channelColor: any;
                detailView: {
                    open: boolean;
                    active: string;
                };
            } & {
                setCurrentAggregation: (samples: number[]) => any;
                setLineUpInput_visibility: (visibility: any) => any;
                setLineUpInput_lineup: (input: any) => any;
                setHoverstate: (state: any, updater: any) => any;
            } & {
                onFilter: any;
            }, "lineUpInput" | "dataset" | "setCurrentAggregation" | "setHoverstate" | "lineUpInput_data" | "lineUpInput_columns" | "currentAggregation" | "channelColor" | "setLineUpInput_lineup" | "setLineUpInput_visibility" | "activeStory" | "pointColorScale" | "detailView">>;
        }[];
    };
};
export declare function CIMEApp(props: CIMEAppProps): JSX.Element;
//# sourceMappingURL=CIMEApp.d.ts.map