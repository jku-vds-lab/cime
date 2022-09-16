/// <reference types="react" />
import { ConnectedProps } from "react-redux";
declare const connector: import("react-redux").InferableComponentEnhancerWithProps<{
    dataset: import("projection-space-explorer").Dataset;
    currentAggregation: {
        aggregation: number[];
        selectedClusters: (string | number)[];
        source: "sample" | "cluster";
    };
    lineUpInput: import("../../State/LineUpInputDuck").LineUpType;
} & {
    setDetailVisibility: (value: any) => any;
    setLineUpInput_filter: (value: any) => any;
}, {}>;
declare type PropsFromRedux = ConnectedProps<typeof connector>;
declare type Props = PropsFromRedux & {
    splitRef: any;
};
export declare const LineUpTabPanel: import("react-redux").ConnectedComponent<({ setDetailVisibility, setLineUpInput_filter, lineUpInput, dataset, currentAggregation, splitRef, }: Props) => JSX.Element, import("react-redux").Omit<{
    dataset: import("projection-space-explorer").Dataset;
    currentAggregation: {
        aggregation: number[];
        selectedClusters: (string | number)[];
        source: "sample" | "cluster";
    };
    lineUpInput: import("../../State/LineUpInputDuck").LineUpType;
} & {
    setDetailVisibility: (value: any) => any;
    setLineUpInput_filter: (value: any) => any;
} & {
    splitRef: any;
}, "lineUpInput" | "dataset" | "currentAggregation" | "setDetailVisibility" | "setLineUpInput_filter">>;
export {};
//# sourceMappingURL=LineUpTabPanel.d.ts.map