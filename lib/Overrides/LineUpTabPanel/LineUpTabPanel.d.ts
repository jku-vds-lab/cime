import { ConnectedProps } from "react-redux";
declare const connector: import("react-redux").InferableComponentEnhancerWithProps<{
    dataset: import("projection-space-explorer").Dataset;
    currentAggregation: {
        aggregation: number[];
        selectedClusters: string[];
        source: "sample" | "cluster";
    };
    lineUpInput: import("../../State/LineUpInputDuck").LineUpType;
} & {
    setLineUpInput_visibility: (value: any) => any;
    setLineUpInput_filter: (value: any) => any;
}, {}>;
declare type PropsFromRedux = ConnectedProps<typeof connector>;
declare type Props = PropsFromRedux & {
    splitRef: any;
};
export declare const LineUpTabPanel: import("react-redux").ConnectedComponent<({ setLineUpInput_visibility, setLineUpInput_filter, lineUpInput, dataset, currentAggregation, splitRef, }: Props) => JSX.Element, import("react-redux").Omit<{
    dataset: import("projection-space-explorer").Dataset;
    currentAggregation: {
        aggregation: number[];
        selectedClusters: string[];
        source: "sample" | "cluster";
    };
    lineUpInput: import("../../State/LineUpInputDuck").LineUpType;
} & {
    setLineUpInput_visibility: (value: any) => any;
    setLineUpInput_filter: (value: any) => any;
} & {
    splitRef: any;
}, "lineUpInput" | "dataset" | "currentAggregation" | "setLineUpInput_visibility" | "setLineUpInput_filter">>;
export {};
//# sourceMappingURL=LineUpTabPanel.d.ts.map